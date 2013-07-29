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
package edu.msu.cme.rdp.graph.filter;

import edu.msu.cme.rdp.kmer.Kmer;

/**
 *
 * @author fishjord
 */
public interface CodonWalker {

    public void jumpTo(char[] s);
    public void jumpTo(Kmer kmer, long fwdHash, long rcHash);

    public NextCodon getNextCodon();
    /**
     * @return alternate amino acid for current codon position.
     * If fails, removes current codon (previous codon becomes current).
     * If it backs up to the starting k-mer, it leaves it in that state
     * and returns 0
     *
     */
    public NextCodon getSibCodon();
    public boolean hasMoreCodons();

    public Byte getNextNucl();
    public Byte getSibNucl();
    public boolean hasMoreNucl();

    public long getFwdHash();
    public long getRcHash();

    /**
     *
     * @return the path starting from the char right after the original kmer,
     *  all the way to the current char pointed by pathPtr
     */
    public String getPathString();

    public int getLength();
}

