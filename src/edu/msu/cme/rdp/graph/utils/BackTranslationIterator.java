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
package edu.msu.cme.rdp.graph.utils;

import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.readseq.utils.NuclBinMapping;
import edu.msu.cme.rdp.readseq.utils.ProteinUtils;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import edu.msu.cme.rdp.seqfilter.filters.ValidAlphabetFilter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * This may say Iterator, but don't expect it to behave exactly like a java
 * iterator
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class BackTranslationIterator {

    private final int[] state;
    private final char[] protMer;
    private final byte[] buf;
    private final int k;
    private final int protLength;
    private boolean hasMore = true;
    private long generated = 1;
    //protected static final Map<Character, List<char[]>> backTranslation;
    protected static final byte[][][] backTrans = new byte[127][][];

    static {
        Map<Character, List<byte[]>> backTranslationMapping = new HashMap();
        final char[] intToNucl = NuclBinMapping.intToChar;
        ProteinUtils.AminoAcid[][][] translTable = ProteinUtils.getInstance().getTranslationTable(11);
        ProteinUtils.AminoAcid aa;
        char p;
        for (byte b1 = 0;b1 < 4;b1++) {
            for (byte b2 = 0;b2 < 4;b2++) {
                for (byte b3 = 0;b3 < 4;b3++) {
                    aa = translTable[b1][b2][b3];
                    p = aa.getAminoAcid();
                    if (!backTranslationMapping.containsKey(p)) {
                        backTranslationMapping.put(Character.toLowerCase(p), new ArrayList());
                        backTranslationMapping.put(Character.toUpperCase(p), new ArrayList());
                    }

                    byte[] add = new byte[]{b1, b2, b3};
                    backTranslationMapping.get(Character.toLowerCase(p)).add(add);
                    backTranslationMapping.get(Character.toUpperCase(p)).add(add);
                }
            }
        }

        for (char residue : SeqUtils.proteinAlphabet) {
            if (!Character.isLowerCase(residue)) {
                continue;
            }

            if(SeqUtils.proteinAmbiguity.contains(residue)) {
                continue;
            }

            if(!backTranslationMapping.containsKey(residue)) {
                //if(residue == 'x') { //ambiguity residue
                System.err.println(residue);
                continue;
            }

            List<byte[]> codonList = backTranslationMapping.get(residue);
            byte[][] codons = new byte[codonList.size()][3];
            for (int index = 0; index < codons.length; index++) {
                codons[index] = codonList.get(index);
            }

            backTrans[residue] = codons;
            backTrans[Character.toUpperCase(residue)] = codons;
        }
    }

    public static long countCombinations(char[] protKmer) {
        long ret = backTrans[protKmer[0]].length;
        for (int index = 1; index < protKmer.length; index++) {
            ret *= backTrans[protKmer[index]].length;
        }

        return ret;
    }

    public BackTranslationIterator(char[] protKmer) {
        k = protKmer.length * 3;
        protLength = protKmer.length;
        state = new int[k];
        buf = new byte[k];
        protMer = protKmer;
        fillBuf();
    }

    public long getGenerated() {
        return generated;
    }

    public byte[] next() {
        if (!hasMore) {
            return null;
        }

        int index = 0;
        for (; index < protLength; index++) {
            state[index]++;
            if (state[index] >= backTrans[protMer[index]].length) {
                state[index] = 0;
                copy(buf, backTrans[protMer[index]][state[index]], index);
            } else {
                copy(buf, backTrans[protMer[index]][state[index]], index);
                break;
            }
        }

        if (index == protMer.length) {
            hasMore = false;
            return null;
        }

        generated++;


        return buf;
    }

    private void fillBuf() {
        for (int index = 0; index < protLength; index++) {
            copy(buf, backTrans[protMer[index]][state[index]], index);
        }
    }

    private static void copy(byte[] dest, byte[] codon, int protIndex) {
        int offset = protIndex * 3;
        dest[offset] = codon[0];
        dest[offset + 1] = codon[1];
        dest[offset + 2] = codon[2];
    }
}
