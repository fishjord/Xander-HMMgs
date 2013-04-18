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
package edu.msu.cme.rdp.graph.cli;

import edu.msu.cme.rdp.graph.utils.BloomFilterStats;
import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.readseq.SequenceFormat;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;

/**
 *
 * @author wangqion
 */
public class BloomFilterBuilder {

    public static class BloomSize {

        public final int bloomSizeLog2;
        public final double fpr;

        public BloomSize(int bloomSizeLog2, double fpr) {
            this.bloomSizeLog2 = bloomSizeLog2;
            this.fpr = fpr;
        }
    }

    public static BloomSize predictBloomSize(File readFile, int kmerSize, int numHashes, int startingM, double desiredFpr) throws IOException {
        long readFileSize = readFile.length();
        long readUntil = (long) (readFileSize * .2);

        SequenceReader reader = new SequenceReader(readFile);
        Sequence seq;
        BloomFilter filter = new BloomFilter(28, numHashes, kmerSize, 8);
        BloomFilter.GraphBuilder builder = filter.new GraphBuilder();

        int seqCount = 0;
        while ((seq = reader.readNextSequence()) != null && reader.getPosition() < readUntil && seqCount < 500000) {
            seqCount++;
            builder.addString(seq.getSeqString().toCharArray());
        }

        int uniqueKmers = (int) filter.getUniqueKmers();
        int predictedKmers = (int) (uniqueKmers * (readFileSize / (double) reader.getPosition()));

        int m = startingM - 1;
        double fpr = Double.MAX_VALUE;

        while (fpr > desiredFpr) {
            m++;
            fpr = Math.pow((1 - Math.pow(Math.E, -numHashes * ((predictedKmers + .5) / (Math.pow(2, m) - 1)))), numHashes);
        }

        return new BloomSize(m, fpr);
    }

    public static void main(String[] args) throws Exception {
        List<File> readFiles = new ArrayList();

        for (int index = 0; index < args.length; index++) {
	    File f = new File(args[index]);
	    if(!f.exists()) {
		break;
	    }
            SequenceFormat format = SeqUtils.guessFileFormat(f);
            if (format == SequenceFormat.UNKNOWN || format == SequenceFormat.EMPTY) {
		break;
            }

            readFiles.add(new File(args[index]));
        }

        args = Arrays.copyOfRange(args, readFiles.size(), args.length);

        if (args.length < 3 || args.length > 5) {
            System.err.println("USAGE: BloomFilterBuilder <read_file> <bloom_out> <kmerSize> <bloomSizeLog2> [# hashCount = 4] [bitsetSizeLog2]");
	    System.err.println("Unexpected number of arguments: " + args.length);
	    System.err.println("Input files: " + readFiles);
            System.exit(1);
        }

        File outputFile = new File(args[0]);

        final int kmerSize = Integer.parseInt(args[1]);
        final int hashSizeLog2 = Integer.parseInt(args[2]);
        final int hashCount;
        final int bitsetSizeLog2;

        if (args.length > 3) {
            hashCount = Integer.parseInt(args[3]);
        } else {
            hashCount = 4;
        }

        if (args.length > 4) {
            bitsetSizeLog2 = Integer.parseInt(args[4]);
        } else {
            if (hashSizeLog2 > 30) {
                bitsetSizeLog2 = 30;
            } else {
                int tmpBitsetSize = 1;
                while (tmpBitsetSize < hashSizeLog2) {
                    tmpBitsetSize <<= 1;
                }

                bitsetSizeLog2 = tmpBitsetSize >> 1;
            }
        }

        if (outputFile.exists()) {
            System.err.println("WARNING: Bloom filter " + outputFile + " already exists, press CTRL+^C to cancel");
        } else {
            outputFile.createNewFile();
        }

        if (!outputFile.canWrite()) {
            throw new IOException("Cannot write to bloom filter file " + outputFile);
        }

        BloomFilter filter = new BloomFilter(hashSizeLog2, hashCount, kmerSize, bitsetSizeLog2);
        BloomFilter.GraphBuilder graphBuilder = filter.new GraphBuilder();

        long seqCount = 0;

        System.err.println("Starting to build bloom filter at " + new Date());
        System.err.println("*  reads file(s):       " + readFiles);
        System.err.println("*  bloom output:     " + outputFile);
        System.err.println("*  kmer size:        " + kmerSize);
        System.err.println("*  hash size log2:   " + hashSizeLog2);
        System.err.println("*  hash count:       " + hashCount);
        System.err.println("*  bitset size log2: " + bitsetSizeLog2);

        long startTime = System.currentTimeMillis();

        for (File readFile : readFiles) {
            SequenceReader reader = new SequenceReader(readFile);
            Sequence seq = null;

            while ((seq = reader.readNextSequence()) != null) {

                seqCount++;
                if ((seqCount % 1000000) == 0) {
                    System.err.println("p: " + seqCount + " kmers added " + graphBuilder.getKmerAdded());
                }

                graphBuilder.addString(seq.getSeqString().toCharArray());
            }
            reader.close();
        }
        long endTime = System.currentTimeMillis();

        ObjectOutputStream oos = new ObjectOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));

        oos.writeObject(filter);
        oos.close();

        BloomFilterStats.printStats(filter, System.out);
        System.err.println("time to build BloomFilter: " + (endTime - startTime) / 60000.0 + " minutes");
    }
}
