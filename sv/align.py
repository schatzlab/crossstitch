import sys
import os
import networkx as nx

reads_file = sys.argv[1]
overlap_file = sys.argv[2]
out_dir = sys.argv[3]

containment_tolerance = 10
permitted_error_pct = 20
    
class read_node_0:  #experimenting on keeping the data in a list rather than a dictionary, not used

    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.hits = [] 
        self.hit_map = None
        
    def add_hit(self, aln_info):
        self.hits.append ( aln_info )

    def __getitem__(self, target_name):
        if self.hit_map == None:
            self.hit_map = dict(zip( [h[0] for h in self.hits], self.hits ) )
        return self.hit_map[target_name]

class read_node:

    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.hit_map = {}
        
    def add_hit(self, aln_info):
        self.hit_map[aln_info[0]] = aln_info

    def __getitem__(self, target_name):
        return self.hit_map[target_name]
    
    @property
    def hits(self):  #another convenient representation of the same data
        return self.hit_map.values()


def get_overlap_data(m4_filename):
    print m4_filename
    overlap_data = {}
    contained_reads = set()
    with open(m4_filename) as m4_f:
        for l in m4_f:
            l = l.strip().split()
            q_name, t_name =l[0:2]
            if q_name == t_name:
                continue
            aln_score = int(l[2])
            aln_idt = float(l[3])
            if aln_idt < 100-permitted_error_pct:
                print 'con'
                continue
            q_strand, q_start, q_end, q_len = ( int(x) for x in l[4:8])
            t_strand, t_start, t_end, t_len = ( int(x) for x in l[8:12])
            if q_len - (q_end - q_start) < containment_tolerance:
                contained_reads.add(q_name)
    
            if t_len - (t_end - t_start) < containment_tolerance:
                contained_reads.add(t_name)
            if q_name not in overlap_data:
                overlap_data[ q_name ] = read_node(q_name, q_len)
            assert q_strand == 0
            if t_name not in [x[0] for x in overlap_data[ q_name ].hits]:
                overlap_data[ q_name ].add_hit( (t_name, 
                                                 aln_score, 
                                                 aln_idt, 
                                                 (q_strand, q_start, q_end, q_len),
                                                 (t_strand, t_start, t_end, t_len) ) )
    
            #symmetrized the alignment record
            if t_name not in overlap_data:
                overlap_data[ t_name ] = read_node(t_name, t_len)
    
            if q_name not in [x[0] for x in overlap_data[ t_name ].hits]:
                if t_strand == 1: 
                    t_start, t_end = t_len - t_end, t_len - t_start
                    q_start, q_end = q_len - q_end, q_len - q_start
                    overlap_data[ t_name ].add_hit( (q_name, 
                                                 aln_score, 
                                                 aln_idt, 
                                                 (0, t_start, t_end, t_len),
                                                 (1, q_start, q_end, q_len) ) )
                else:
                    overlap_data[ t_name ].add_hit( (q_name, 
                                                 aln_score, 
                                                 aln_idt, 
                                                 (0, t_start, t_end, t_len),
                                                 (0, q_start, q_end, q_len) ) )
    return overlap_data, contained_reads

def generate_overlap_gml(overlap_data, contained_reads, gml_filename):
    G=nx.DiGraph()
    node_in_graph = set()
    print contained_reads
    for q_name in overlap_data:
        if q_name in contained_reads:
            continue
        if q_name not in node_in_graph:
            G.add_node(q_name)
        targets = overlap_data[ q_name ].hits
        print targets
        targets_3prime = [ h for h in targets if h[4][1] < containment_tolerance and h[0] not in contained_reads]
        targets_5prime = [ h for h in targets if h[3][1] < containment_tolerance and h[0] not in contained_reads]
        print str(targets_3prime) + ' ' + str(targets_5prime)
        targets_3prime.sort(key = lambda k:k[1])
        targets_5prime.sort(key = lambda k:k[1])

        if len(targets_3prime) > 0:
            t = targets_3prime[0]
            t_name = t[0]
            if t_name not in node_in_graph:
                G.add_node(t_name)
            G.add_edge(q_name, t_name)

        if len(targets_5prime) > 0:
            t = targets_5prime[0]
            t_name = t[0]
            if t_name not in node_in_graph:
                G.add_node(t_name)
            G.add_edge(q_name, t_name)
    nx.write_gml(G, gml_filename)
    
def filter_overlap_hit(overlap_hit):
    containment_tolerance = 50
    q_strand, q_start, q_end, q_len = overlap_hit[3]
    t_strand, t_start, t_end, t_len = overlap_hit[4]
    if q_len - (q_end - q_start) < containment_tolerance and  t_len - (t_end - t_start) / t_len < containment_tolerance:
        return False
    elif overlap_hit[0] in contained_reads:
        return False
    else:
        return True

def get_best_overlap_graph(overlap_data):
    best_overlap_graph = {}
    containment_tolerance = 50
    
    for q_name in overlap_data:
        if q_name in contained_reads:
            continue
        
        if q_name not in best_overlap_graph:
            best_overlap_graph[q_name] = {"5p":None, "3p":None}
        
        
        targets = overlap_data[ q_name ].hits
        
        targets_5prime = [ h for h in targets if h[3][1] < containment_tolerance and filter_overlap_hit(h)]
        targets_3prime = [ h for h in targets if h[4][1] < containment_tolerance and filter_overlap_hit(h)]
    
        targets_5prime.sort(key = lambda k:k[1])
        targets_3prime.sort(key = lambda k:k[1])
    
        if len(targets_5prime) > 0:
            best_overlap_graph[q_name]["5p"] = targets_5prime[0]
            t_name = targets_5prime[0][0]
    
                
        if len(targets_3prime) > 0:
            best_overlap_graph[q_name]["3p"] = targets_3prime[0]
            t_name = targets_3prime[0][0]
    
    return best_overlap_graph
    
def proper_overlap_hit(overlap_hit):
    containment_tolerance = 50
    q_strand, q_start, q_end, q_len = overlap_hit[3]
    t_strand, t_start, t_end, t_len = overlap_hit[4]
    if q_start < containment_tolerance:
        if t_len - t_end > containment_tolerance:
            return False
        else:
            return True
    if t_start < containment_tolerance:
        if q_len - q_end > containment_tolerance:
            return False
        else:
            return True
        
def find_path( q_name_end, frag_used = set() ):
    reverse_end = {"3p":"5p", "5p":"3p"}
    path = []
    path_q_name = set()
    q_name, end = q_name_end 
    if end == "5p":
        reversing_path = True
    else:
        reversing_path = False
        
    while 1:
        if q_name in frag_used:
            if reversing_path:
                path.reverse()
            return path, "frag_used"
        
        path.append( (q_name, end) )
        path_q_name.add(q_name)
        if q_name not in best_overlap_graph:
            if reversing_path:
                path.reverse()
            return path, "terminate_1" 
        next_hit = best_overlap_graph[q_name][end]
        #print next_hit
        if next_hit == None:
            if reversing_path:
                path.reverse()
            return path, "terminate_2"
        
        if next_hit[0] in best_overlap_graph: #if not mutual good hits, break the path
            
            # Using mutual best hit might be to strigent, 
            # bh = []
            #if best_overlap_graph[next_hit[0]]["5p"]:
            #    bh.append( best_overlap_graph[next_hit[0]]["5p"][0] )
            #if best_overlap_graph[next_hit[0]]["3p"]:
            #    bh.append( best_overlap_graph[next_hit[0]]["3p"][0] )
            
            bh = [h[0] for h in overlap_data[next_hit[0]].hits if proper_overlap_hit(h)]
       
            if q_name not in bh:
                if reversing_path:
                    path.reverse()
                return path, "branch"
        
        q_name = next_hit[0]
        
        if q_name in path_q_name:
            if reversing_path:
                path.reverse()
            return path, "circle"
        
        if next_hit[4][0] == 1: #reversed strand
            end = reverse_end[end]
            
def rc_seq(seq):
    rev_map = dict(zip("acgtACGTNn-","tgcaTGCANn-"))
    return "".join([rev_map[c] for c in seq[::-1]])

def layout_path(full_path, frag_used, out_fn, tig_name):
    
    if len(full_path) == 0 or full_path[0][0] in frag_used:
        return None
    
    if len(full_path) == 1:
        with open(out_dir + "/" + "singleton_"+out_fn, "w") as out_seq_file:
            seq = seq_db[full_path[0][0]]
            print >>out_seq_file, ">%s" % ( "singleton_" + tig_name  )
            print >>out_seq_file, seq
            frag_used.add(full_path[0][0])
        return None
    
    first = full_path[0]
    offset = 0
    current_orientation = first[1]
    revserse_orientation = {"5p":"3p", "3p":"5p"}
    tig_seq = []
    frag_in_layout = set()
    frag_in_layout.add(first[0])
    for second in full_path[1:]:
        overlap_hit = overlap_data[first[0]][second[0]]
        
        q_strand, q_start, q_end, q_len = overlap_hit[3]
        t_strand, t_start, t_end, t_len = overlap_hit[4]
        #print overlap_hit, offset, current_orientation
        seq = seq_db[first[0]]
        if current_orientation == "3p":
            seq = rc_seq(seq)
        del tig_seq[offset:-1]
        tig_seq.extend(seq)
        
        if current_orientation == "5p":
            offset += q_start
        elif current_orientation == "3p":
            offset += q_len - q_end
        if t_strand == 1:
           current_orientation  = revserse_orientation[current_orientation] 
        frag_in_layout.add(first[0])
        first = second
        if second[0] in frag_in_layout:
            break
        if second[0] in frag_used:
            break

    seq = seq_db[first[0]]
    if current_orientation == "3p":
        seq = rc_seq(seq)
    del tig_seq[offset:-1]
    tig_seq.extend(seq)
    frag_in_layout.add(first[0])
    
    frag_used.update( frag_in_layout )
    
    tig_seq = tig_seq[0:offset+len(seq)]
    tig_seq = "".join(tig_seq)
    with open(out_dir + "/" + out_fn, "w") as out_seq_file:
        print >>out_seq_file, ">%s" % tig_name
        print >>out_seq_file, tig_seq
    
overlap_data, contained_reads = get_overlap_data(overlap_file)
best_overlap_graph = get_best_overlap_graph(overlap_data)

seq_db = {}
count = 0
with open(reads_file) as seq_file:
    for l in seq_file:
        l = l.strip()
        if len(l) == 0:
            continue
        if ' ' in l:
            l = l[0:l.index(' ')]
        if l[0] == ">":
            count += 1
            name = l[1:]
            continue
        else:
            seq_db[name] = l

                
frag_used = set()

all_best_overlap_frags = set(best_overlap_graph.keys())

unused_frag = all_best_overlap_frags - frag_used
i = 0

if len(unused_frag) == 0:
    for q_name in seq_db:
        out_fn = "tig_%05d.fa" % i
        tig_name ="tig_%05d" % i
        with open(out_dir + "/" + "singleton_"+out_fn, "w") as out_seq_file:
            seq = seq_db[q_name]
            print >>out_seq_file, ">%s" % ( "singleton_" + tig_name  )
            print >>out_seq_file, seq
        i += 1


while len(unused_frag) != 0:

    len_qname = [ (overlap_data[x].length, x) for x in unused_frag  ]
    len_qname.sort()
    q_name = len_qname[-1][1]
    print "iteration: ",i
    print "frag is used?", q_name in frag_used
    if q_name in frag_used:
        continue

    path_5p, s_5p = find_path( (q_name, "5p"), frag_used  )
    path_3p, s_3p = find_path( (q_name, "3p"), frag_used  )
    print "seeding frag:", q_name, "Length:", overlap_data[q_name].length
    print "number of unused frag:", len(unused_frag), "total overlapped frag:", len(path_3p)+len(path_5p)
    print len(path_5p), s_5p, len(path_3p), s_3p
    print "--------------------------"
    #print path_5p[0], path_3p[-1]
    if len(path_5p) + len(path_3p) == 0:
        out_fn = "tig_%05d.fa" % i
        tig_name ="tig_%05d" % i
        with open(out_dir + "/" + "singleton_"+out_fn, "w") as out_seq_file:
            seq = seq_db[q_name]
            print >>out_seq_file, ">%s" % ( "singleton_" + tig_name  )
            print >>out_seq_file, seq
        frag_used.add(q_name)
        unused_frag = all_best_overlap_frags - frag_used
        i += 1
        continue
        
    if len(path_5p) > 0:
        #assert path_5p[-1] == path_3p[0]

        full_path = path_5p + path_3p[1:]
    else:
        full_path = path_3p
    layout_path(full_path, frag_used, "tig_%05d.fa" % i, "tig_%05d" % i)
    unused_frag = all_best_overlap_frags - frag_used
    i += 1

