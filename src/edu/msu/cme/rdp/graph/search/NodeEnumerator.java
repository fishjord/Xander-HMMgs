
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

import edu.msu.cme.rdp.alignment.hmm.MostProbableHCostHMM;
import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.graph.filter.NextCodon;
import static edu.msu.cme.rdp.alignment.hmm.TSC.*;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import edu.msu.cme.rdp.graph.filter.PathHolder;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.readseq.SequenceType;
import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class NodeEnumerator {

    private double hweight = 1;
    private AStarNode next;
    private NextCodon protEmission = null;
    private Byte nextNucl = null;
    private char emission;
    private double matchTrans;
    private double insTrans;
    private double delTrans;
    private Kmer nextKmer;
    private int nextState;
    private boolean newKmer;
    private final ProfileHMM hmm;
    private final boolean protSearch;
    private final MostProbableHCostHMM hcost;
    private final long mask = (1L << 61) - 1;

    public NodeEnumerator(ProfileHMM hmm) {
        this.hmm = hmm;
        this.protSearch = hmm.getAlphabet() == SequenceType.Protein;
        this.hcost = hmm.getHCost();
    }

    public Set<AStarNode> enumerateNodes(AStarNode curr, CodonWalker walker, Set<Kmer> seenKmers) {
        Set<AStarNode> ret = new HashSet();

        nextState = curr.stateNo + 1;

        switch (curr.state) {
            case 'm':
                matchTrans = hmm.tsc(curr.stateNo, MM);
                insTrans = hmm.tsc(curr.stateNo, MI);
                delTrans = hmm.tsc(curr.stateNo, MD);
                break;
            case 'd':
                matchTrans = hmm.tsc(curr.stateNo, DM);
                insTrans = Double.NEGATIVE_INFINITY;
                delTrans = hmm.tsc(curr.stateNo, DD);
                break;
            case 'i':
                matchTrans = hmm.tsc(curr.stateNo, IM);
                insTrans = hmm.tsc(curr.stateNo, II);
                delTrans = Double.NEGATIVE_INFINITY;
                break;
            default:
                throw new RuntimeException("I hate you.");
        }


        walker.jumpTo(curr.kmer, curr.fwdHash, curr.rcHash);

        if (protSearch) {
            protEmission = walker.getNextCodon();
        } else {
            nextNucl = walker.getNextNucl();
        }

        double maxMatchEmission = hmm.getMaxMatchEmission(nextState);
        while ((protSearch && protEmission != null) || (!protSearch && nextNucl != null)) {
            if (protSearch) {
                if (protEmission.getAminoAcid() == '*') {
                    protEmission = walker.getSibCodon();
                    continue;
                }

                int codon = protEmission.getCodon() & 127;
                byte b1 = (byte) (codon & 0x3);
                byte b2 = (byte) (codon >> 2 & 0x3);
                byte b3 = (byte) (codon >> 4 & 0x3);
                nextKmer = curr.kmer.shiftLeft(b3).shiftLeft(b2).shiftLeft(b1);
                emission = protEmission.getAminoAcid();
            } else {
                nextKmer = curr.kmer.shiftLeft((byte) (nextNucl & 3));
                emission = NuclBinMapping.intToChar[nextNucl];
            }

            final long fwdHash = walker.getFwdHash();
            final long rcHash = walker.getRcHash();

            /**
             * ************************************
             *
             * MATCH NODE
             *
             *************************************
             */
            newKmer = !seenKmers.contains(nextKmer);
            next = new AStarNode(curr, nextKmer, fwdHash, rcHash, nextState, 'm');

            next.realScore = curr.realScore + matchTrans + hmm.msc(nextState, emission);
            next.emission = emission;
            next.thisNodeScore = matchTrans + hmm.msc(nextState, emission) - maxMatchEmission;
            next.length = curr.length + 1;
            next.score = (curr.score + next.thisNodeScore);
            next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + hweight * hcost.computeHeuristicCost('m', nextState)));
            next.hasNewKmer = curr.hasNewKmer || newKmer;
            next.indels = curr.indels;

            ret.add(next);

            /**
             * ************************************
             *
             * INSERT NODE
             *
             *************************************
             */
            if (curr.state != 'd') { //Transitions from delete to insert aren't allowed, don't waste time computing stuff
                next = new AStarNode(curr, nextKmer, fwdHash, rcHash, curr.stateNo /*
                         * Inserts don't advance the state
                         */, 'i');
                next.realScore = curr.realScore + insTrans + hmm.isc(nextState, emission);
                next.emission = emission;
                next.thisNodeScore = insTrans + hmm.isc(nextState, emission);
                next.length = curr.length + 1;
                next.score = (curr.score + next.thisNodeScore);
                next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + hweight * hcost.computeHeuristicCost('i', curr.stateNo)));
                next.hasNewKmer = curr.hasNewKmer || newKmer;
                next.indels = curr.indels + 1;

                ret.add(next);
            }

            if (protSearch) {
                protEmission = walker.getSibCodon();
            } else {
                nextNucl = walker.getSibNucl();
            }
        }

        /**
         * ************************************
         *
         * DELETE NODE
         *
         *************************************
         */
        if (curr.state != 'i') {
            next = new AStarNode(curr, curr.kmer, curr.fwdHash, curr.rcHash, nextState, 'd');

            next.realScore = curr.realScore + delTrans;
            next.emission = '-';
            next.thisNodeScore = delTrans - maxMatchEmission;
            next.length = curr.length;
            next.score = (curr.score + next.thisNodeScore);
            next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + hweight * hcost.computeHeuristicCost('d', nextState)));
            next.hasNewKmer = curr.hasNewKmer;
            next.indels = curr.indels + 1;

            ret.add(next);
        }

        return ret;
    }
}
