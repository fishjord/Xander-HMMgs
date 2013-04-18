package edu.msu.cme.rdp.graph.search;

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


import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.readseq.SequenceType;

/**
 *
 * @author fishjord
 */
public class SearchTarget {

    private String kmer;
    private int frame;
    private int startState;
    private ProfileHMM forwardHmm;
    private ProfileHMM reverseHmm;
    private BloomFilter filter;
    private boolean prot;
    private String geneName, querySeqid, refSeqid;

    public SearchTarget(String geneName, String querySeqid, String refSeqid, String kmer, int frame, int startState, ProfileHMM forwardHmm, ProfileHMM reverseHmm, BloomFilter filter) {
        this.kmer = kmer;
        this.frame = frame;
        this.startState = startState;
        this.forwardHmm = forwardHmm;
        this.reverseHmm = reverseHmm;
        this.filter = filter;
        this.prot = forwardHmm.getAlphabet() == SequenceType.Protein;
        this.geneName = geneName;
        this.querySeqid = querySeqid;
        this.refSeqid = refSeqid;
    }

    public String getGeneName() {
        return geneName;
    }

    public String getQuerySeqid() {
        return querySeqid;
    }

    public String getRefSeqid() {
        return refSeqid;
    }

    public int getFrame() {
        return frame;
    }

    public ProfileHMM getForwardHmm() {
        return forwardHmm;
    }

    public ProfileHMM getReverseHmm() {
        return reverseHmm;
    }

    public String getKmer() {
        return kmer;
    }

    public int getStartState() {
        return startState;
    }

    public BloomFilter getFilter() {
        return filter;
    }

    public boolean isProt() {
        return prot;
    }
}
