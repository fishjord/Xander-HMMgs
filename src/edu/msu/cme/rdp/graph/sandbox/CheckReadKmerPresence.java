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
package edu.msu.cme.rdp.graph.sandbox;

import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import edu.msu.cme.rdp.kmer.trie.KmerGenerator;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import java.io.*;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author fishjord
 */
public class CheckReadKmerPresence {

    public static void printStats(BloomFilter filter, PrintStream out) {


        long n = filter.getUniqueKmers();
        long m = (long) Math.pow(2, filter.getHashSizeLog2());
        int k = filter.getHashCount();

        //(1-e^(-k*((n+.5)/(m-1))))^k
        double falsePositiveRate = Math.pow((1 - Math.pow(Math.E, -k * ((n + .5) / (m - 1)))), k);

        out.println("Bloom filter created on:       " + filter.getCreatedOn());
        out.println("Serializable id:               " + BloomFilter.serialVersionUID);
        out.println();
        out.println("Bloom filter size log 2:       " + filter.getHashSizeLog2());
        out.println("Bloom filter size (bits) (m):  " + m);
        out.println();
        out.println("Number of bitsets:             " + filter.getNumBitsets());
        out.println("Bitset size (bits):            " + filter.getBitsetSize());
        out.println("Bitset size log2:              " + filter.getBitsetSizeLog2());
        out.println();
        out.println("Number of hashes (k):          " + filter.getHashCount());
        out.println("Hash function name:            " + filter.getHasherClassName());
        out.println();
        out.println("Bitset Mask:                   " + StringUtils.leftPad(Long.toBinaryString(filter.getBitsetMask()), 64, '0'));
        out.println("Hash mask:                     " + StringUtils.leftPad(Long.toBinaryString(filter.getHashMask()), 64, '0'));
        out.println();
        out.println("Kmer length:                   " + filter.getKmerSize());
        out.println();
        out.println("Total kmers in bloom filter:   " + filter.getTotalKmers());
        out.println("Total strings inserted:        " + filter.getTotalStrings());
        out.println("Total unique kmers:            " + filter.getUniqueKmers());
        out.println("Predicted false positive rate: " + falsePositiveRate);
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            System.err.println("USAGE: CheckReadKmerPresence <bloom_filter> <nucl_seq_file>");
            System.exit(1);
        }

        File bloomFile = new File(args[0]);
        SeqReader reader = new SequenceReader(new File(args[1]));

        ObjectInputStream ois = ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(bloomFile)));
        BloomFilter filter = (BloomFilter) ois.readObject();
        ois.close();

        printStats(filter, System.out);
        Sequence seq;
        CodonWalker walker = null;
        while ((seq = reader.readNextSequence()) != null) {
            int kmerNum = 0;
            for (char[] kmer : KmerGenerator.getKmers(seq.getSeqString(), filter.getKmerSize())) {
                System.out.print(seq.getSeqName() + "\t" + (++kmerNum) + "\t" + kmer + "\t");
                try {
                    walker = filter.new RightCodonFacade(kmer);
                    System.out.println("true");
                } catch (Exception e) {
                    System.out.println("false");
                }

            }
        }
    }
}
