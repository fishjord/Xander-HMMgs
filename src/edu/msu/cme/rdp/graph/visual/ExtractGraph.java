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
package edu.msu.cme.rdp.graph.visual;

import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import edu.msu.cme.rdp.graph.filter.InvalidDNABaseException;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.set.NuclKmerGenerator;
import edu.msu.cme.rdp.kmer.trie.KmerGenerator;
import edu.msu.cme.rdp.readseq.SequenceFormat;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.utils.SeqUtils;
import java.io.*;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class ExtractGraph {

    private final BloomFilter bloom;
    private final int radius;
    private final Graph graph;
    private int kmersProcessed = 0;

    public ExtractGraph(BloomFilter bloom, int radius) {
        this.bloom = bloom;
        this.radius = radius;
        graph = new Graph();
    }

    /**
     * I'm deliberately processing kmers multiple times if they show up multiple
     * times because it's hard to know if we've already processed this kmer at
     * this depth, one option to speed this up would be to keep track of the
     * depth and kmer and not process duplicates, we'll see how slow this is
     *
     * @param kmer
     */
    public void processKmer(String kmerStr) {
        processKmer(kmerStr.toCharArray());
    }

    public void processKmer(char[] kmerStr) {
        Kmer kmer = new Kmer(kmerStr);
        CodonWalker walker;

        walker = bloom.new LeftCodonFacade(kmerStr);
        processKmer(walker, kmer, radius, false);

        walker = bloom.new RightCodonFacade(kmerStr);
        processKmer(walker, kmer, radius, true);

        kmersProcessed++;
    }

    private void processKmer(CodonWalker walker, Kmer kmer, int depth, boolean fwd) {
        if (depth <= 0) {
            return;
        }

        Byte b = walker.getNextNucl();
        if (b == null) {
            return;
        }

        Kmer next;
        do {

            if (fwd) {
                next = kmer.shiftLeft(b);
            } else {
                next = kmer.shiftRight(b);
            }

            graph.connect(kmer, next);

            processKmer(walker, next, depth - 1, fwd);

        } while ((b = walker.getSibNucl()) != null);
    }

    public Graph getGraph() {
        return graph;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 3 && args.length != 4) {
            System.err.println("USAGE: ExtractGraph <bloom_filter> <kmer_or_seqfile> <radius> [diam]");
            System.exit(1);
        }

        File bloomInFile = new File(args[0]);
        File inputKmerFile = new File(args[1]);
        int radius = Integer.valueOf(args[2]);
        int diam = 0;

        if(args.length == 4) {
            diam = Integer.valueOf(args[2]);
        }
        long startTime = System.currentTimeMillis();

        BloomFilter filter;
        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(bloomInFile)));
        filter = (BloomFilter) ois.readObject();
        ois.close();
        System.err.println("Bloomfilter loaded in " + (System.currentTimeMillis() - startTime) / 1000.0f + "s");

        int k = filter.getKmerSize();
        ExtractGraph graphExtractor = new ExtractGraph(filter, radius);
        startTime = System.currentTimeMillis();

        SequenceFormat format = SeqUtils.guessFileFormat(inputKmerFile);
        if (format == SequenceFormat.UNKNOWN) {
            BufferedReader reader = new BufferedReader(new FileReader(inputKmerFile));
            String kmer;

            while ((kmer = reader.readLine()) != null) {
                if (kmer.length() != k) {
                    System.err.println(kmer + " length doesn't match bloomfilter k=" + k);
                    continue;
                }
                graphExtractor.processKmer(kmer);
            }
        } else {
            SeqReader reader = new SequenceReader(inputKmerFile);
            Sequence seq;
            NuclKmerGenerator kmerGen;

            while((seq = reader.readNextSequence()) != null) {

                for(char[] kmer : KmerGenerator.getKmers(seq.getSeqString(), k)) {
                    try {
                        graphExtractor.processKmer(kmer);
                    } catch(InvalidDNABaseException e) {
                    }
                }
            }

        }

	System.err.println("Graph loaded from " + args[0] + " in " + (System.currentTimeMillis() - startTime) / 1000.0f + "s");
	System.err.println("Vertices: " + graphExtractor.getGraph().getNumVertices());
	System.err.println("Edges: " + graphExtractor.getGraph().getNumEdges());

	int componentCount = 1;
	for(Graph comp : graphExtractor.getGraph().getComponents(diam)) {
	    comp.writeDot(new File("component_" + componentCount + ".gv"));
	    System.err.println("\tWriting out component " + componentCount + " vertex count=" + comp.getNumVertices() + " edge count=" + comp.getNumEdges());
	    componentCount++;
	}
    }
}
