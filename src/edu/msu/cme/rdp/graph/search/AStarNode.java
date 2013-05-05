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

import edu.msu.cme.rdp.graph.filter.PathHolder;
import edu.msu.cme.rdp.kmer.Kmer;
import java.io.Serializable;

/**
 *
 * @author fishjord
 */
public class AStarNode implements Serializable, Comparable<AStarNode> {

    public boolean partial = true;
    public AStarNode discoveredFrom;
    public final Kmer kmer;
    public final long fwdHash, rcHash;
    public double score;
    public final char state;
    public final int stateNo;
    public int fval;
    public boolean hasNewKmer;
    public double thisNodeScore;
    public int indels;
    public double realScore;
    public int length;
    public char emission;

    public AStarNode(AStarNode discoveredFrom, Kmer kmer, long fwdHash, long rcHash, int stateNo, char state) {
        this.discoveredFrom = discoveredFrom;
        this.fwdHash = fwdHash;
        this.rcHash = rcHash;
        this.kmer = kmer;
        this.stateNo = stateNo;
        this.state = state;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final AStarNode other = (AStarNode) obj;

        if (!this.kmer.equals(other.kmer)) {
            return false;
        }
        if (this.state != other.state) {
            return false;
        }
        if (this.stateNo != other.stateNo) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 5;

        hash = 13 * kmer.hashCode();
        hash = 13 * hash + this.state;
        hash = 13 * hash + this.stateNo;
        return hash;
    }

    @Override
    public String toString() {
        return ((state == 'd') ? "-" : kmer.toString()) + " " + state + " " + stateNo + " " + " " + score + " " + fval;
    }

    public int compareTo(AStarNode o) {

        int ret = o.fval - fval;

        if (ret == 0) {
            ret = stateNo - o.stateNo;
        }

        if (ret == 0) {
            char[] k1 = kmer.toString().toCharArray();
            char[] k2 = o.kmer.toString().toCharArray();

            for (int index = 0; index < k1.length; index++) {
                if (k1[index] != k2[index]) {
                    ret = k2[index] - k1[index];
                    break;
                }
            }

        }

        if (ret == 0) {
            int s1 = 0;
            int s2 = 0;
            switch (state) {
                case 'm':
                    s1 = 3;
                    break;
                case 'd':
                    s1 = 2;
                    break;
                case 'i':
                    s1 = 1;
                    break;
            }
            switch (o.state) {
                case 'm':
                    s2 = 3;
                    break;
                case 'd':
                    s2 = 2;
                    break;
                case 'i':
                    s2 = 1;
                    break;
            }

            ret = s2 - s1;
        }

        return ret;
    }
}
