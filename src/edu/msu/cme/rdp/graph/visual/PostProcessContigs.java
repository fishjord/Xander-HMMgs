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
package edu.msu.cme.rdp.graph.visual;

import edu.msu.cme.rdp.graph.search.AStarNode;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.trie.KmerGenerator;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author fishjord
 */
public class PostProcessContigs {

    public static void main(String[] args) throws IOException {
        if(args.length != 2 && args.length != 3) {
            System.err.println("USAGE: PostProcessContigs <nucl_contigs> <k> [min_diameter=0]");
            System.exit(1);
        }

        SeqReader reader = new SequenceReader(new File(args[0]));
        int k = Integer.valueOf(args[1]);
        int diam = 0;
        if(args.length == 3) {
            diam = Integer.valueOf(args[2]);
        }

        Sequence seq;
        Graph graph = new Graph();
	long startTime = System.currentTimeMillis();

        while((seq = reader.readNextSequence()) != null) {
            AStarNode node = null;
            for(char[] kmer : KmerGenerator.getKmers(seq.getSeqString(), k)) {
                AStarNode next = new AStarNode(node, new Kmer(kmer), 0, 0, 0, 'm');
                node = next;
            }

            graph.connectAll(node);
        }

        reader.close();

	System.err.println("Graph loaded from " + args[0] + " in " + (System.currentTimeMillis() - startTime) / 1000.0f + "s");
	System.err.println("Vertices: " + graph.getNumVertices());
	System.err.println("Edges: " + graph.getNumEdges());


	int componentCount = 1;
	for(Graph comp : graph.getComponents(diam)) {
	    comp.writeDot(new File("component_" + componentCount + ".gv"));
	    System.err.println("\tWriting out component " + componentCount + " vertex count=" + comp.getNumVertices() + " edge count=" + comp.getNumEdges());
	    componentCount++;
	}
    }
}
