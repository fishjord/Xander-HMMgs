package edu.msu.cme.rdp.graph.search;

import edu.msu.cme.rdp.kmer.io.KmerStart;

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
/**
 *
 * @author fishjord
 */
public class SearchResult {
    public static enum SearchDirection { left, right };

    private String kmer;
    private String nuclSeq;
    private String protSeq;
    private long time;
    private double nats;
    private double bits;
    private SearchDirection searchDirection;
    private int mpos;
    private String alignSeq;
    private SearchTarget start;

    public SearchResult(SearchTarget start, String kmer, String nuclSeq, String alignSeq, String protSeq, SearchDirection dir, int mpos, double nats, double bits, long time) {
        this.bits = bits;
        this.nats = nats;
        this.time = time;
        this.mpos = mpos;

        this.kmer = kmer;
        this.nuclSeq = nuclSeq;
        this.protSeq = protSeq;
        this.searchDirection = dir;
        this.alignSeq = alignSeq;
        this.start = start;
    }

    public SearchTarget getStart() {
        return start;
    }

    public String getKmer() {
        return kmer;
    }

    public String getNuclSeq() {
        return nuclSeq;
    }

    public String getProtSeq() {
        return protSeq;
    }

    public long getTime() {
        return time;
    }

    public double getBits() {
        return bits;
    }

    public double getNats() {
        return nats;
    }

    public SearchDirection getSearchDirection() {
        return searchDirection;
    }

    public int getModelPosition() {
        return mpos;
    }

    public String getAlignSeq() {
        return alignSeq;
    }
}
