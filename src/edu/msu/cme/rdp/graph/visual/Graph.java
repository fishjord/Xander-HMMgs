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
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;

/**
 *
 * @author fishjord
 */
public class Graph implements Serializable {
    protected Map<Kmer, Set<Kmer>> edges = new HashMap();
    protected boolean combined;

    public void connect(Kmer v1, Kmer v2) {
	if(!edges.containsKey(v1)) {
	    edges.put(v1, new HashSet());
	}

	if(!edges.containsKey(v2)) {
	    edges.put(v2, new HashSet());
	}

	edges.get(v1).add(v2);
	edges.get(v2).add(v1);
    }

    public int getNumVertices() {
	return edges.size();
    }

    public int getNumEdges() {
	int ret = 0;
	for(Set<Kmer> s : edges.values()) {
	    ret += s.size();
	}

	return ret;
    }

    public void connectAll(AStarNode node) {
	while(node != null && node.discoveredFrom != null) {
	    connect(node.kmer, node.discoveredFrom.kmer);
	    node = node.discoveredFrom;
	}
    }

    public List<Graph> getComponents(int minDiam) {
        List<Graph> ret = new ArrayList();
        if (edges.isEmpty()) {
            return ret;
        }

        Set<Kmer> allNodes = new HashSet(edges.keySet());

        while (!allNodes.isEmpty()) {
            Kmer next = allNodes.iterator().next(); //Grab any node, don't care
            Graph comp = new Graph();
	    Set<Kmer> thisComp = new HashSet();

            int d = connectedTo(next, thisComp, comp);

	    if(d > minDiam) {
		ret.add(comp);
	    }
            allNodes.removeAll(thisComp);
        }

        return ret;
    }

    private int connectedTo(Kmer kmer, Set<Kmer> thisComp, Graph comp) {
        if (thisComp.contains(kmer)) {
            return 0;
        }

	thisComp.add(kmer);
	int ret = 0;

        for (Kmer next : edges.get(kmer)) {
	    comp.connect(kmer, next);
            ret = Math.max(ret, connectedTo(next, thisComp, comp));
        }

	return ret + 1;
    }

    public void writeDot(File outFile) throws IOException {
	PrintWriter out = new PrintWriter(outFile);
	Map<Kmer, Integer> labels = new HashMap();
	try {
	    out.println("graph sg {");
	    for (Kmer ls : edges.keySet()) {
		int label = labels.size();
		labels.put(ls, label);
		out.println("\t" + ls + "[label=\"" + label + "\"];");
		for (Kmer rs : edges.get(ls)) {
		    if(labels.containsKey(rs)) {
			continue;
		    }

		    out.println("\t" + ls + " -- " + rs + ";");
		}
		out.println();
	    }
	    out.println("}");

	} finally {
	    out.close();
	}
    }
}
