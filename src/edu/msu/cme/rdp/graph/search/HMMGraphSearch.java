/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.graph.search;

import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.alignment.hmm.scoring.ForwardScorer;
import edu.msu.cme.rdp.alignment.hmm.scoring.HMMScorer;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import edu.msu.cme.rdp.graph.filter.PathHolder;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import java.io.IOException;
import java.util.*;

/**
 *
 * @author fishjord
 */
public class HMMGraphSearch {

    public static class PartialResult {

        public String maxSeq;
        public String alignment;
        public double maxScore;
        public long searchTime;
    }

    public static class HackTerminateException extends RuntimeException {
    }
    public static final int INT_SCALE = 10000; //This is the number of sigfigs in a HMMER3 model, so it works out quite nicely if you ask me
    private static final int upperBound = Integer.MIN_VALUE;
    private final int maxk;
    //private PrintStream openedKmerStream;
    //private PrintStream closedKmerStream;

    public HMMGraphSearch(int maxk) {
        this.maxk = maxk;
        /*
         * try { this.openedKmerStream = new
         * PrintStream("all_opened_kmers.txt"); this.closedKmerStream = new
         * PrintStream("all_closed_kmers.txt"); } catch(IOException ignored) {
         *
         * }
         */
    }

    public List<SearchResult> search(SearchTarget target) throws InterruptedException {
        String framedKmer = target.getKmer();
        int frame = target.getFrame();
        List<SearchResult> ret = new ArrayList();

        int lStartingState = target.getReverseHmm().M() - target.getStartState() - target.getKmer().length() / 3 + 2;

        List<PartialResult> leftParts = kpathsSearch(target.getReverseHmm(), lStartingState, framedKmer, target.getFilter().new LeftCodonFacade(target.getKmer()), false);
        for (PartialResult r : leftParts) {
            String nuclSeq = r.maxSeq + framedKmer;
            String alignment = r.alignment + framedKmer.toUpperCase();
            String protSeq = null;
            String scoringSeq = nuclSeq;

            HMMScorer scorer = new ForwardScorer(target.getForwardHmm(), -1);

            if (target.isProt()) {
                scoringSeq = protSeq = ProteinUtils.getInstance().translateToProtein(nuclSeq, true, 11);
            }

            for (char c : scoringSeq.toCharArray()) {
                scorer.consume(c);
            }

            ret.add(new SearchResult(target, target.getKmer(), nuclSeq, alignment, protSeq, SearchResult.SearchDirection.left, lStartingState, r.maxScore, scorer.getMaxScore(), r.searchTime));
        }

        List<PartialResult> rightParts = kpathsSearch(target.getForwardHmm(), target.getStartState(), framedKmer, target.getFilter().new RightCodonFacade(target.getKmer()), true);

        for (PartialResult r : rightParts) {
            String nuclSeq = framedKmer + r.maxSeq;
            String alignment = framedKmer.toUpperCase() + r.alignment;
            String protSeq = null;
            String scoringSeq = nuclSeq;

            HMMScorer scorer = new ForwardScorer(target.getForwardHmm(), -1);

            if (target.isProt()) {
                scoringSeq = protSeq = ProteinUtils.getInstance().translateToProtein(nuclSeq, true, 11);
            }

            for (char c : scoringSeq.toCharArray()) {
                scorer.consume(c);
            }

            ret.add(new SearchResult(target, target.getKmer(), nuclSeq, alignment, protSeq, SearchResult.SearchDirection.right, target.getStartState(), r.maxScore, scorer.getMaxScore(), r.searchTime));
        }

        return ret;
    }

    public List<AStarNode> searchGraph(SearchTarget target) throws InterruptedException {
        String framedKmer = target.getKmer();
        int frame = target.getFrame();
        List<AStarNode> ret = new ArrayList();

        int lStartingState = target.getReverseHmm().M() - target.getStartState() - target.getKmer().length() / 3 + 2;

        List<CandidatePath> leftParts = kpathsSearchGraph(target.getReverseHmm(), lStartingState, framedKmer, target.getFilter().new LeftCodonFacade(target.getKmer()), false);

        for (CandidatePath r : leftParts) {
            ret.add(r.get(r.length() - 1));
        }

        List<CandidatePath> rightParts = kpathsSearchGraph(target.getForwardHmm(), target.getStartState(), framedKmer, target.getFilter().new RightCodonFacade(target.getKmer()), true);

        for (CandidatePath r : rightParts) {
            ret.add(r.get(r.length() - 1));
        }

        return ret;
    }

    private List<PartialResult> kpathsSearch(ProfileHMM hmm, int startingState, String framedWord, CodonWalker walker, boolean forward) throws InterruptedException {
        List<PartialResult> results = new ArrayList();
        for (CandidatePath path : kpathsSearchGraph(hmm, startingState, framedWord, walker, forward)) {
            AStarNode goal = path.get(path.length() - 1);
            results.add(partialResultFromGoal(goal, forward, hmm.getAlphabet() == SequenceType.Protein, framedWord.length(), path.generationTime));
        }

        return results;
    }

    private List<CandidatePath> kpathsSearchGraph(ProfileHMM hmm, int startingState, String framedWord, CodonWalker walker, boolean forward) throws InterruptedException {
        List<CandidatePath> bestPaths = new ArrayList();
        PriorityQueue<CandidatePath> candidatePaths = new PriorityQueue<CandidatePath>();
        Map<AStarNode, Set<AStarNode>> shortestPathEdges = new HashMap();
        Set<Kmer> seenKmers = new HashSet();


        long kTime = System.currentTimeMillis();
        //PrintStream out = new PrintStream(forward? "right.txt" : "left.txt");
        try {
            AStarNode goalNode = astarSearch(hmm, startingState, framedWord, walker, forward, seenKmers, new HashSet());

            CandidatePath bestPath = new CandidatePath(goalNode, seenKmers);
            bestPath.generationTime = (System.currentTimeMillis() - kTime);
            bestPaths.add(bestPath);
            //out.println(bestPath.score + "\t" + bestPath.k + "\t" + bestPath.i + "\t" + partialResultFromGoal(goalNode, forward, hmm.getAlphabet() == SequenceType.Protein, 0).maxSeq);

            while (bestPaths.size() < maxk) {   //Where k is the current kth shortest path
                CandidatePath pathAk = bestPaths.get(bestPaths.size() - 1);
                kTime = System.currentTimeMillis();

                /*
                 * Candidate generation
                 */
                for (int i = pathAk.i; i < pathAk.length() - 1; i++) { //Lawler's Observation
                    CandidatePath base = pathAk.subpath(i + 1);

                    AStarNode starting = base.get(i);
                    AStarNode ak_i_1 = pathAk.get(i + 1);

                    if (shortestPathEdges.containsKey(starting) && shortestPathEdges.get(starting).contains(ak_i_1)) {
                        break;
                    }

                    if (!shortestPathEdges.containsKey(starting)) {
                        shortestPathEdges.put(starting, new HashSet());
                    }

                    shortestPathEdges.get(starting).add(ak_i_1);
                    goalNode = astarSearch(hmm, starting, walker, seenKmers, shortestPathEdges.get(starting));

                    int before = seenKmers.size();
                    CandidatePath spur = new CandidatePath(goalNode, seenKmers);

                    CandidatePath candidate = spur;
                    candidate.i = i;
                    candidate.k = bestPaths.size();

                    if (!bestPaths.contains(candidate) && seenKmers.size() > before) {
                        candidatePaths.add(candidate);
                    } else {
                        shortestPathEdges.get(starting).remove(ak_i_1);
                    }

                    //out.println(candidate.score + "\t" + candidate.k + "\t" + candidate.i + "\t" + partialResultFromGoal(goalNode, forward, hmm.getAlphabet() == SequenceType.Protein, 0).maxSeq);
                }

                CandidatePath kthPath = candidatePaths.poll();

                if (kthPath == null || Double.isInfinite(kthPath.score)) {
                    break;
                }
                kthPath.generationTime = (System.currentTimeMillis() - kTime);

                bestPaths.add(kthPath);
            }
        } catch (HackTerminateException e) {
            System.err.println("Terminated on path " + bestPaths.size() + " (candidates=" + candidatePaths.size() + ")");
            //openedKmerStream.close();
            //closedKmerStream.close();

            throw e;
        } catch (IOException ignore) {
        }

        return bestPaths;
    }

    /**
     *
     * NOTE: Walker -MUST- be initialized to the passed starting kmer
     *
     * @param hmm
     * @param startingState
     * @param framedWord
     * @param walker
     * @param forward
     * @param seenKmers
     * @param disallowedLinks
     * @return
     * @throws IOException
     */
    private AStarNode astarSearch(final ProfileHMM hmm,
            int startingState, String framedWord,
            CodonWalker walker,
            boolean forward,
            Set<Kmer> seenKmers,
            Set<AStarNode> disallowedLinks) throws IOException, InterruptedException {
        framedWord = framedWord.toLowerCase();

        char[] startingCodon = framedWord.substring(framedWord.length() - 3).toCharArray();
        Kmer kmer;
        if (forward) {
            kmer = new Kmer(framedWord.toCharArray());
        } else {
            kmer = new Kmer(new StringBuilder(framedWord).reverse().toString().toCharArray());
        }

        AStarNode startingNode;
        if (hmm.getAlphabet() == SequenceType.Protein) {
            startingNode = new AStarNode(null, kmer, walker.getFwdHash(), walker.getRcHash(), startingState + (framedWord.length() / 3) - 1, 'm');
        } else {
            startingNode = new AStarNode(null, kmer, walker.getFwdHash(), walker.getRcHash(), startingState - 1, 'm');
        }

        startingNode.fval = 0;
        startingNode.score = 0;

        return astarSearch(hmm, startingNode, walker, seenKmers, disallowedLinks);
    }

    private AStarNode astarSearch(final ProfileHMM hmm,
            AStarNode startingNode,
            CodonWalker walker,
            Set<Kmer> seenKmers,
            Set<AStarNode> disallowedLinks) throws IOException, InterruptedException {
        NodeEnumerator nodeEnumerator = new NodeEnumerator(hmm);
        PriorityQueue<AStarNode> open = new PriorityQueue<AStarNode>();
        Set<AStarNode> closed = new HashSet();
        AStarNode curr;

        int nodesClosed = 0;
        int nodesOpened = 0;
        int window = 25;

        double[] bestPerBaseNats = new double[hmm.M() - startingNode.stateNo];
        Arrays.fill(bestPerBaseNats, Double.NEGATIVE_INFINITY);
        double perBaseNats = 0;

        Thread currThread = Thread.currentThread();
        long start = System.currentTimeMillis();
        int maxIndels = 5;//(int)(hmm.M() * .05 + .5);

        if (startingNode.stateNo >= hmm.M()) {   //Huh...well I guess we don't get much choice in the matter now do we?
            return startingNode;
        }
        //PrintStream out = new PrintStream((forward ? "forward" : "reverse") + ".gv");
        //out.println("digraph " + (forward ? "forward" : "reverse") + "{");

        //First step, enumerate all the nodes and remove any disallowed transitions
        //This way we only have to look at the set (disallowedLinks) once instead of
        //during every iteration (which was silly)
        for (AStarNode next : nodeEnumerator.enumerateNodes(startingNode, walker, seenKmers)) {
            if (!disallowedLinks.contains(next)
                    && next.score > upperBound) {
                nodesOpened++;
                open.add(next);
            }
        }
        //Decide the intermediate goal
        AStarNode interGoal = (open.peek() == null) ? startingNode : open.peek();
        int max = 10000;
        AStarNode o;
        double mem;
        long t;

        //While we have more things to close
        while ((curr = open.poll()) != null) {
            if (closed.contains(curr)) {
                continue;
            }

            nodesClosed++;

            if (curr.stateNo >= hmm.M()) { //We're at an "end" state
                if (curr.hasNewKmer) {  //If it has a new kmer, great
                    return curr;
                } else { //Otherwise move on
                    continue;
                }
            }

            /*
             * if (interGoal.score > curr.score || interGoal.stateNo <
             * curr.stateNo) { //If we've got a better interim result keep track
             * of it interGoal = curr; }
             */

            closed.add(curr);

            if ((closed.size() % 250000) == 0) {
                mem = getMemRatio();
                //System.err.println("Open set size: " + open.size() + ", Closed: " + closed.size() + ", node: " + curr + " mem ratio: " + mem + " time: " + (System.currentTimeMillis() - start) / 1000.0f);
                if (mem > .75) {
                    t = System.currentTimeMillis();
                    open.removeAll(closed);
                    System.gc();
                    //    System.err.println("\tMemory reclaimation time: " + (System.currentTimeMillis() - t) / 1000.0f + "s, mem ratio after reclaimation: " + getMemRatio());
                }

                if(Thread.interrupted()) {
                    throw new InterruptedException();
                }
            }

            /*
             * int mlen = curr.stateNo - startingNode.stateNo; if (mlen > 0) {
             * perBaseNats = curr.score / mlen;//getAverageScore(curr, window);
             * if (perBaseNats > bestPerBaseNats[mlen]) { bestPerBaseNats[mlen]
             * = perBaseNats; } else if (perBaseNats < bestPerBaseNats[mlen] *
             * 1.05) { //System.out.println("Evicting " + curr + " because it's
             * per base score (" + perBaseNats + ") from " +
             * startingNode.stateNo + " to " + curr.stateNo + " (" + mlen + ")
             * is not better than 85% of the highest score (" +
             * bestPerBaseNats[mlen] + ")"); continue; } }
             */

            //Look at the adjacent nodes
            for (AStarNode next : nodeEnumerator.enumerateNodes(curr, walker, seenKmers)) {
                //Make sure we haven't already seen something better
                if (next.score > upperBound && next.indels < maxIndels) {
                    nodesOpened++;
                    open.add(next);
                }
            }
        }
        interGoal.score = Double.NEGATIVE_INFINITY;
        interGoal.fval = Integer.MIN_VALUE;
        return interGoal;
    }

    /**
     * Assemble the path, from the goal to the start
     */
    public static PartialResult partialResultFromGoal(AStarNode goal, boolean forward, boolean protSearch, int kmerLength, long searchTime) {
        StringBuilder nuclSeq = new StringBuilder();
        StringBuilder alignmentSeq = new StringBuilder();

        char[] gap = (protSearch) ? new char[]{'-', '-', '-'} : new char[]{'-'};

        PartialResult result = new PartialResult();
        result.maxScore = goal.score;

        if (goal != null) {
            while (goal.discoveredFrom != null) {
                char[] kmer = goal.kmer.toString().toCharArray();
                char[] emission;

                if (protSearch) {
                    if(forward) {
                        //If we're to the right the codon is at the end of the kmer
                        emission = new char[]{kmer[kmer.length - 3], kmer[kmer.length - 2], kmer[kmer.length - 1]};
                    } else {
                        //If we're going to the left the codon is still at the right
                        //end of this kmer BUT it is in reverse order (DIFFERENT
                        //than in the path returned by CodonWalker.getPathString())
                        emission = new char[]{kmer[kmer.length - 1], kmer[kmer.length - 2], kmer[kmer.length - 3]};
                    }
                } else {
                    //In the single emission case the last character in the kmer
                    //is always right
		    emission = new char[]{kmer[kmer.length - 1]};
                }

                if (goal.state == 'd') {
                    if (forward) {
                        alignmentSeq.insert(0, gap);
                    } else {
                        alignmentSeq.append(gap);
                    }
                } else if (goal.state == 'm') {
                    if (forward) {
                        alignmentSeq.insert(0, new String(emission).toUpperCase());
                    } else {
                        alignmentSeq.append(new String(emission).toUpperCase());
                    }
                } else if (goal.state == 'i') {
                    if (forward) {
                        alignmentSeq.insert(0, new String(emission).toLowerCase());
                    } else {
                        alignmentSeq.append(new String(emission).toLowerCase());
                    }
                }


                if (goal.state != 'd') { //No emission on delete states
                    if (forward) {
                        //prepend for forward
                        nuclSeq.insert(0, emission);
                    } else {
                        //append for reverse (we're building in the 'right' direction)
                        nuclSeq.append(emission);
                    }
                }

                goal = goal.discoveredFrom;
            }
        }
        result.maxSeq = nuclSeq.toString();
        result.alignment = alignmentSeq.toString();
        result.searchTime = searchTime;

        return result;
    }

    private static double getAverageScore(AStarNode node, int over) {
        double ret = 0;
        int index;
        for (index = 0; index < over && node != null; index++, node = node.discoveredFrom) {
            ret += node.thisNodeScore;
        }

        if (index != over) {
            throw new IllegalArgumentException("Path does not go back " + over + " nodes");
        }

        return ret / over;
    }

    public static double getMemRatio() {
        return (double) (Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / Runtime.getRuntime().maxMemory();
    }
}
