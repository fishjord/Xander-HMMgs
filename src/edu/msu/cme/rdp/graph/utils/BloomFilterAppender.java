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
package edu.msu.cme.rdp.graph.utils;

import edu.msu.cme.rdp.graph.utils.BloomFilterStats;
import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.Date;

/**
 *
 * @author wangqion
 */
public class BloomFilterAppender {

    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            System.err.println("USAGE: BloomFilterAppender <bloomfilter> <read_file>");
            System.exit(1);
        }

        File bloomFilter = new File(args[0]);
        BloomFilter filter = BloomFilter.fromFile(bloomFilter);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();

        long seqCount = 0;

        args = Arrays.copyOfRange(args, 1, args.length);
        System.err.println("Starting to build bloom filter at " + new Date());
        System.err.println("*  reads file(s):       " + Arrays.asList(args));
        System.err.println("*  bloom output:     " + bloomFilter);
        System.err.println("*  kmer size:        " + filter.getKmerSize());
        System.err.println("*  hash size log2:   " + filter.getHashSizeLog2());
        System.err.println("*  hash count:       " + filter.getHashCount());
        System.err.println("*  bitset size log2: " + filter.getBitsetSize());

        long startTime = System.currentTimeMillis();

        for (String f : args) {
            File readFile = new File(f);
            SequenceReader reader = new SequenceReader(readFile);
            Sequence seq;

            while ((seq = reader.readNextSequence()) != null) {

                seqCount++;
                if ((seqCount % 1000000) == 0) {
                    System.err.println("p: " + seqCount + " kmers added " + graphBuilder.getKmerAdded());
                }

                graphBuilder.addString(seq.getSeqString().toCharArray());
            }
            reader.close();
        }

        BloomFilterStats.printStats(filter, System.out);
        long endTime = System.currentTimeMillis();

        System.err.println("time to build BloomFilter: " + (endTime - startTime) / 60000.0 + " minutes");

        ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(bloomFilter.getAbsolutePath() + ".appended")));

        oos.writeObject(filter);
        oos.close();
    }
}
