
/*
 * Copyright (C) 2013 Jordan Fish <fishjord at msu.edu>
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
package edu.msu.cme.rdp.graph.sandbox;

import com.sun.tools.internal.xjc.util.NullStream;
import edu.msu.cme.rdp.alignment.hmm.HMMER3bParser;
import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.alignment.hmm.TSC;
import edu.msu.cme.rdp.alignment.hmm.scoring.HMMScorer;
import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import edu.msu.cme.rdp.graph.search.AStarNode;
import edu.msu.cme.rdp.graph.search.NodeEnumerator;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class ExploreRenderCombinedGraph {

    private final BloomFilter bloom;
    private final ProfileHMM forHmm;
    private final ProfileHMM revHmm;
    private final int maxRadius;
    private final boolean allowGaps;
    private final boolean prot;
    private boolean forward;

    public ExploreRenderCombinedGraph(BloomFilter bloom, ProfileHMM forHmm, ProfileHMM revHmm, int radius, boolean allowGaps) {
        this.bloom = bloom;
        this.forHmm = forHmm;
        this.revHmm = revHmm;
        this.maxRadius = radius;
        this.allowGaps = allowGaps;
        this.prot = forHmm.getAlphabet() == SequenceType.Protein;
    }

    public void explore(String kmer, int hmmState, PrintStream out) {
        Set<AStarNode> exploredNodes = new HashSet();
        Map<AStarNode, Set<AStarNode>> edges = new HashMap();
        out.println("digraph " + kmer + " {");

        forward = true;
        CodonWalker walker = bloom.new RightCodonFacade(kmer);
        AStarNode startingNode = getStartNode(kmer, hmmState, forHmm, walker);
        out.println("\t" + getNodeLabel(startingNode, 0));
        walk(walker, startingNode, out, new NodeEnumerator(forHmm), forHmm);

        forward = false;
        walker = bloom.new LeftCodonFacade(kmer);
        startingNode = getStartNode(kmer, forHmm.M() - hmmState - kmer.length() / ((prot)? 3 : 1), revHmm, walker);
        walk(walker, startingNode, out, new NodeEnumerator(revHmm), revHmm);

        out.println("}");
    }

    private void walk(CodonWalker walker, AStarNode start, PrintStream out, NodeEnumerator enumerator, ProfileHMM hmm) {
        Stack<AStarNode> toVisit = new Stack();
        Set<AStarNode> seenNodes = new HashSet();

        toVisit.push(start);

        System.err.println("#radius\tnum_residues\tstate_no\tstate\temission\tscore\tbits");
        double best = 0;
        while (!toVisit.isEmpty()) {
            AStarNode prev = toVisit.pop();
	    if(prev.stateNo >= hmm.M()) {
		break;
	    }
            int radius = prev.indels;
            char[] kmer = prev.kmer.toString().toCharArray();
            //System.err.println(new String(kmer));

            walker.jumpTo(prev.kmer, prev.fwdHash, prev.rcHash);
            char[] emission;

            if (prot) {
                if (forward) {
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
            double score = prev.realScore + Math.log(2.0 / (prev.length + 2)) * 2;
            if(score > -50) {
            //    best = score;
                System.err.println(radius + "\t" + prev.length + "\t" + prev.stateNo + "\t" + prev.state + "\t" + new String(emission) + "\t" + ProteinUtils.getInstance().translateToProtein(new String(emission), true, 11) + "\t" + prev.score + "\t" + prev.realScore + "\t" + score + "\t" + (score - HMMScorer.getNull1(prev.length)) / HMMScorer.ln2);
            }


            for (AStarNode next : enumerator.enumerateNodes(prev, walker, new HashSet())) {
                /*if(next.state != 'm') {
                 if(!allowGaps || radius < 20) {
                 continue;
                 }
                 }*/

                if (!allowGaps && next.state != 'm') {
                    continue;
                }

                //if ((next.realScore - HMMScorer.getNull1(next.length)) / HMMScorer.ln2 < -10) {
                //    continue;
                //}

                if (radius < maxRadius) {
		    if(!seenNodes.contains(next)) {
			seenNodes.add(next);
			toVisit.push(next);
			next.indels = radius + 1;
			out.println(getNodeLabel(next, radius + 1));
		    }
		    out.println(getEdge(prev, next));
                }
            }
        }
    }

    private String getId(AStarNode node) {
        StringBuilder ret = new StringBuilder();
        char[] kmer = node.kmer.toString().toCharArray();

        ret.append(kmer).append("_").append(node.stateNo).append(node.state);
        return ret.toString();
    }

    private String getEdge(AStarNode node1, AStarNode node2) {
        return getId(node1) + "->" + getId(node2) + ";";
    }

    private String getNodeLabel(AStarNode node, int radius) {
        StringBuilder ret = new StringBuilder();
        char[] kmer = node.kmer.toString().toCharArray();

        ret.append(getId(node));

        ret.append("[style=\"filled\",shape=\"box\",fillcolor=\"");
        if (node.state == 'm') {
            ret.append("#3399FF");
        } else if (node.state == 'd') {
            ret.append("#FF3300");
        } else if (node.state == 'i') {
            ret.append("#ff00FF");
        } else {
            throw new RuntimeException("I hate you.");
        }
        char[] emission;

        if (prot) {
            if (forward) {
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

        ret.append("\",label=\"");
        ret.append(String.format("Emission: %s\\nState: %d\\nRadius: %d\\nScore: %.2f \\nbits: %.2f", (emission.length == 1 ? emission[0] : ProteinUtils.getInstance().translateToProtein(new String(emission), true, 11)), node.stateNo, radius, node.score, (node.realScore - HMMScorer.getNull1(node.length)) / HMMScorer.ln2));
        ret.append("\"];");

        return ret.toString();
    }

    private AStarNode getStartNode(String startingKmer, int startingState, ProfileHMM hmm, CodonWalker walker) {
        Kmer kmer = new Kmer(startingKmer.toCharArray());

        AStarNode startingNode;
        if (prot) {
            startingNode = new AStarNode(null, kmer, walker.getFwdHash(), walker.getRcHash(), startingState + startingKmer.length() / 3, 'm');
        } else {
            startingNode = new AStarNode(null, kmer, walker.getFwdHash(), walker.getRcHash(), startingState, 'm');
        }

        String k = startingKmer;
        if(forward && prot) {
            k = ProteinUtils.getInstance().translateToProtein(startingKmer, false, 11);
        } else if(!forward && prot) {
            StringBuilder tmp = new StringBuilder(ProteinUtils.getInstance().translateToProtein(startingKmer, false, 11));
            k = tmp.reverse().toString();
        }


        startingNode.length = 7;
        startingNode.realScore = realScoreStart(hmm, k, startingState);
        startingNode.fval = 0;
        startingNode.score = scoreStart(hmm, k, startingState);


        return startingNode;
    }

    private float scoreStart(ProfileHMM hmm, String startingKmer, int startingState) {
        System.err.println("Starting kmer: " + startingKmer + ", state: " + startingState);
        float ret = 0;
        char[] residues = startingKmer.toCharArray();
        for (int index = 1; index <= residues.length && index + startingState < hmm.M() - 1; index++) {
            ret += hmm.msc(startingState + index, residues[index - 1]) + hmm.tsc(startingState + index, TSC.MM) - hmm.getMaxMatchEmission(startingState + index);
            System.err.println(index + "\t" + residues[index - 1] + "\t" + ret);
        }
        System.err.println("Starting word score: " + ret);


        return ret;
    }

    private float realScoreStart(ProfileHMM hmm, String startingKmer, int startingState) {
        System.err.println("Starting kmer: " + startingKmer + ", state: " + startingState);
        float ret = 0;
        char[] residues = startingKmer.toCharArray();
        for (int index = 1; index <= residues.length && index + startingState < hmm.M(); index++) {
            ret += hmm.msc(startingState + index, residues[index - 1]) + hmm.tsc(startingState + index -1, TSC.MM);
            System.err.println(index + "\t" + residues[index - 1] + "\t" + ret + "\t" + (ret + Math.log(2.0 / (index + startingState + 2)) * 2 - HMMScorer.getNull1(index + startingState)) / HMMScorer.ln2);
        }
        System.err.println("Starting word score: " + ret);


        return ret;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 6 && args.length != 7) {
            System.err.println("USAGE: ExploreRenderCombinedGraph <bloom> <startingKmer> <hmmState> <forHmmModel> <revHmmModel> <radius> [allowGaps=true]");
            System.exit(1);
        }

        File bloomFile = new File(args[0]);
        String startingKmer = args[1];
        int startingState = Integer.valueOf(args[2]);
        File forHmmFile = new File(args[3]);
        File revHmmFile = new File(args[4]);
        int radius = Integer.valueOf(args[5]);
        boolean allowGaps = true;
        Kmer kmer = new Kmer(startingKmer.toCharArray());

        if (args.length == 7) {
            allowGaps = Boolean.valueOf(args[6]);
        }

        ProfileHMM forHmm = HMMER3bParser.readModel(forHmmFile);
        ProfileHMM revHmm = HMMER3bParser.readModel(revHmmFile);
        boolean prot = forHmm.getAlphabet() == SequenceType.Protein;

        BloomFilter bloom;
        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(bloomFile)));
        bloom = (BloomFilter) ois.readObject();
        ois.close();

        ExploreRenderCombinedGraph explorer = new ExploreRenderCombinedGraph(bloom, forHmm, revHmm, radius, allowGaps);

	PrintStream out = new PrintStream("combined_graph.dot");
        explorer.explore(startingKmer, startingState, out);
	out.close();
    }
}
