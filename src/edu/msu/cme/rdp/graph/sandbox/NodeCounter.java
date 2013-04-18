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
import edu.msu.cme.rdp.kmer.io.KmerStart;
import edu.msu.cme.rdp.kmer.io.KmerStartsReader;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class NodeCounter {

    public static int countNodes(CodonWalker walker, Set<Long> seenHashes, boolean prot, int radius) {
        if(radius <= 0) {
            return 0;
        }

        int ret = 1;

        if(prot) {
            if(walker.getNextCodon() != null) {
                ret += countNodes(walker, seenHashes, prot, radius - 1);
            }

            while(walker.getSibCodon() != null) {
                if(walker.getNextCodon() != null) {
                    ret += countNodes(walker, seenHashes, prot, radius - 1);
                }
            }
        } else {
            if(walker.getNextNucl()!= null) {
                ret += countNodes(walker, seenHashes, prot, radius - 1);
            }

            while(walker.getSibNucl() != null) {
                if(walker.getNextNucl() != null) {
                    ret += countNodes(walker, seenHashes, prot, radius - 1);
                }
            }
        }

        return ret;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 3 && args.length != 4) {
            System.err.println("USAGE: NodeCounter <bloom_filter> <kmer_starts> <model_length> [prot = true]");
            System.exit(1);
        }

        File bloomFile = new File(args[0]);
        File kmerStarts = new File(args[1]);
        int radius = Integer.valueOf(args[2]);
        boolean prot = true;
        if(args.length == 4) {
            prot = Boolean.valueOf(args[3]);
        }

        long startTime;

        startTime = System.currentTimeMillis();
        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(bloomFile)));
        BloomFilter bloom = (BloomFilter) ois.readObject();
        ois.close();
        System.err.println("Bloom filter loaded in " + (System.currentTimeMillis() - startTime) + " ms");

        System.err.println("Starting node counting search at " + new Date());
        System.err.println("*  Kmer file:               " + kmerStarts);
        System.err.println("*  Bloom file:              " + bloomFile);
        System.err.println("*  Searching prot?:         " + prot);
        System.err.println("*  radius:                  " + radius);

        startTime = System.currentTimeMillis();

        int kmerCount = 0;
        int count;
        CodonWalker walker;

        KmerStart line;
        KmerStartsReader reader = new KmerStartsReader(kmerStarts);
        long timer;

        try {
            while ((line = reader.readNext()) != null) {
                kmerCount++;

                if (line.getMpos() == 0) {
                    System.err.println("Skipping line " + line);
                    continue;
                }


                timer = System.currentTimeMillis();
                walker = bloom.new LeftCodonFacade(line.getNuclKmer());
                count = countNodes(walker, new HashSet(), prot, Math.min(line.getMpos(), radius));
                System.out.println(line.getNuclKmer() + "\t" + line.getMpos() + "\tleft\t" + count + "\t" + (System.currentTimeMillis() - timer) / 1000.0 + "s");

                timer = System.currentTimeMillis();
                walker = bloom.new RightCodonFacade(line.getNuclKmer());
                count = countNodes(walker, new HashSet(), prot, radius - line.getMpos());
                System.out.println(line.getNuclKmer()+ "\t" + line.getMpos() + "\tright\t" + count + "\t" + (System.currentTimeMillis() - timer) / 1000.0 + "s");
            }

            System.err.println("Finished processing " + kmerCount + " kmers in " + (System.currentTimeMillis() - startTime) / 1000.0 + "s");
        } finally {
            reader.close();
        }

    }
}
