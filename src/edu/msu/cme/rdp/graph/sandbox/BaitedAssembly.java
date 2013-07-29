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
package edu.msu.cme.rdp.graph.sandbox;

import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import edu.msu.cme.rdp.kmer.trie.KmerGenerator;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import edu.msu.cme.rdp.graph.filter.InvalidDNABaseException;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author Jordan Fish <fishjord at msu.edu>
 */
public class BaitedAssembly {

    private final BloomFilter bloom;
    private final BloomFilter.GraphState state;
    private static final int maxDepth = 1000;

    public BaitedAssembly(BloomFilter bloom) {
        this.bloom = bloom;
        state = bloom.new GraphState();
    }

    private void dfs(CodonWalker walker, List<String> assemblies, Set<Long> visited) {
        Byte base = walker.getNextNucl();
        long fwd, rc, hash;

        if (base == null || walker.getLength() > maxDepth) {
            if(walker.getLength() > 0) {
                assemblies.add(walker.getPathString());
            }
            return;
        }

        fwd = walker.getFwdHash();
        rc = walker.getRcHash();
        hash = (fwd < rc) ? fwd : rc;

        if(visited.contains(hash)) {
            return;
        }

        visited.add(hash);

        do {
            dfs(walker, assemblies, visited);
        } while(walker.getSibNucl() != null);
    }

    public List<String> assemble(Sequence bait) {
        char[] kmer = null;

        //We want to add the bait sequence to the assembly, so we need to know which kmer we started assembling from
        //The assembly algorithm does not return the starting kmer, so the starting offset is the end of the first kmer
        int offset = bloom.getKmerSize();
        int finalOffset = offset;
        for (char[] k : new KmerGenerator(bait.getSeqString(), bloom.getKmerSize())) {

            try {
                state.setState(k);
                if (state.hasCurrent()) {
                    finalOffset = offset;
                    kmer = k;
                }
            } catch (InvalidDNABaseException e) {
            }
            offset++;
        }

        List<String> assemblies = new ArrayList();

        if (kmer == null) {
            return assemblies;
        }

        CodonWalker walker = bloom.new RightCodonFacade(kmer);
        Set<Long> visited = new HashSet();

        dfs(walker, assemblies, visited);

        List<String> ret = new ArrayList();
        for(String assembly : assemblies) {
            ret.add(bait.getSeqString().substring(0, finalOffset) + assembly);
        }

        return ret;
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            System.err.println("USAGE: BaitedAssembly <bloom_filter> <bait_seqs>");
            System.exit(1);
        }

        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(args[0])));
        BloomFilter bloom = (BloomFilter) ois.readObject();
        ois.close();

        SeqReader reader = new SequenceReader(new File(args[1]));
        Sequence seq;

        BaitedAssembly assembler = new BaitedAssembly(bloom);

        FastaWriter out = new FastaWriter(System.out);
        int contigs;
        long startTime;

        while ((seq = reader.readNextSequence()) != null) {
            contigs = 1;
            startTime = System.currentTimeMillis();
            for (String assembly : assembler.assemble(seq)) {
                out.writeSeq(seq.getSeqName() + "_contig_" + contigs, assembly);
                System.err.println(seq.getSeqName() + "\t" + contigs + "\t" + assembly.length() + "\t" + (System.currentTimeMillis() - startTime) / 1000.0f);
                contigs++;
            }
        }

        out.close();

        reader.close();
    }
}
