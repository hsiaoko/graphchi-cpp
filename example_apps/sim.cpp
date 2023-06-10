
/**
 * @file
 * @author  Aapo Kyrola <akyrola@cs.cmu.edu>
 * @version 1.0
 *
 * @section LICENSE
 *
 * Copyright [2012] [Aapo Kyrola, Guy Blelloch, Carlos Guestrin / Carnegie Mellon University]
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 
 *
 * @section DESCRIPTION
 *
 * Application for computing the connected components of a graph.
 * The algorithm is simple: on first iteration each vertex sends its
 * id to neighboring vertices. On subsequent iterations, each vertex chooses
 * the smallest id of its neighbors and broadcasts its (new) label to
 * its neighbors. The algorithm terminates when no vertex changes label.
 *
 * @section REMARKS
 *
 * This application is interesting demonstration of the asyncronous capabilities
 * of GraphChi, improving the convergence considerably. Consider
 * a chain graph 0->1->2->...->n. First, vertex 0 will write its value to its edges,
 * which will be observed by vertex 1 immediatelly, changing its label to 0. Nexgt,
 * vertex 2 changes its value to 0, and so on. This all happens in one iteration.
 * A subtle issue is that as any pair of vertices a<->b share an edge, they will
 * overwrite each others value. However, because they will be never run in parallel
 * (due to deterministic parallellism of graphchi), this does not compromise correctness.
 *
 * @author Aapo Kyrola
 */

#include <stdlib.h>

#include <cmath>
#include <string>

#include "graphchi_basic_includes.hpp"
#include "util/labelanalysis.hpp"
#include "util/bitmap.hpp"

using namespace graphchi;

int         iterationcount = 0;
bool        scheduler = false;

/**
 * Type definitions. Remember to create suitable graph shards using the
 * Sharder-program. 
 */
typedef vid_t VertexDataType;       // vid_t is the vertex id type
typedef vid_t EdgeDataType;

struct Vertex {
    unsigned vid_;
    size_t indegree_ = 0;
    size_t outdegree_ = 0;
    unsigned label_;
    unsigned *out_edges_ = nullptr;
};

struct ImmutableCSR {
    size_t num_vertexes_ = 0;
    size_t num_edges_ = 0;
    Vertex ** vertexes_ = nullptr;
};

void ShowImmutableCSR(ImmutableCSR & graph){
    printf("num_vertexes_: %d, num_edges: %d\n", graph.num_vertexes_, graph.num_edges_);
    printf(" out: ");
    for(size_t i =0; i < graph.num_vertexes_; i++){        
        printf("%d, ", graph.vertexes_[i]->vid_);
    }
    printf("  \n  ");
}

struct MatchSet {
      Bitmap* indicator_ = nullptr;
      Bitmap** sim_sets_ = nullptr;
      size_t x_ = 0;
      size_t y_ = 0;
            
      MatchSet(const size_t x, const size_t y, const bool init = false) {
          x_ = x;
          y_ = y;
          indicator_ = new Bitmap(x);
          indicator_->clear();
          sim_sets_ = (Bitmap**)malloc(sizeof(Bitmap*) * x);
          
          if (!init)
              for (size_t i = 0; i < x; i++) sim_sets_[i] = nullptr;
          else {
              for (size_t i = 0; i < x; i++) {
                  sim_sets_[i] = new Bitmap(y);
                  sim_sets_[i]->clear();
              }
          }
      }

      ~MatchSet() {
          if (indicator_ != nullptr) delete indicator_;
          if (sim_sets_ != nullptr) {
              for (size_t i = 0; i < x_; i++)
                  if (sim_sets_[i] == nullptr) delete sim_sets_[i];
              free(sim_sets_);
          }
      }
};

ImmutableCSR* InitPatternNClique(size_t n) {
    ImmutableCSR * graph = new ImmutableCSR;
    graph->num_vertexes_ = n;
    graph->num_edges_ = n * (n - 1);
    graph->vertexes_ = (Vertex **)malloc(sizeof(Vertex*) * graph->num_vertexes_);
    for(size_t i = 0; i < graph->num_vertexes_; i++) {
       Vertex *vertex = new Vertex;
       vertex->vid_ = i;
       vertex->label_ = rand() % 5;
       vertex->indegree_ = n-1;
       vertex->outdegree_ = n-1;
       vertex->out_edges_ = (unsigned *)malloc(sizeof(unsigned) * (n - 1));
       for(size_t j = 0; j < graph->num_vertexes_; j++){
           if(i == j) continue;
           printf("vid: %d, add edge %d", i, j);
           vertex->out_edges_[j] = i;
       }
       printf("\n");
       graph->vertexes_[i] = vertex;
    }
    return graph;
}


/**
 * GraphChi programs need to subclass GraphChiProgram<vertex-type, edge-type> 
 * class. The main logic is usually in the update function.
 */
struct SimPragram : public GraphChiProgram<VertexDataType, EdgeDataType> {
    
    bool converged;

    ImmutableCSR * pattern_;
    MatchSet * match_set_;
    /**
     *  Vertex update function.
     *  On first iteration ,each vertex chooses a label = the vertex id.
     *  On subsequent iterations, each vertex chooses the minimum of the neighbor's
     *  label (and itself). 
     */

    bool kernel_check_childs(graphchi_vertex<VertexDataType, EdgeDataType> & u, unsigned pattern_vid){
        
        size_t count = 0;
        for(int i=0; i < u.num_edges(); i++) {
            vid_t nbr_u_label = u.edge(i)->get_data();
            auto nbr_u_id = u.edge(i)->vertex_id();

            for(size_t j = 0; j < pattern_->vertexes_[pattern_vid]->outdegree_; j++)
            {
                auto nbr_v_id= pattern_->vertexes_[j]->vid_;
                auto nbr_v_label= pattern_->vertexes_[j]->label_;

                if(nbr_u_label == nbr_v_label)
                {
                    if(match_set_->indicator_->get_bit(nbr_v_id) == 0) continue;
                    count++;
                }
            
            }
        
            return count== pattern_->vertexes_[pattern_vid]->outdegree_;
    
        }
    }

    bool CheckU(graphchi_vertex<VertexDataType, EdgeDataType> & u, unsigned pattern_vid) {
        auto uid = u.vertexid;
        return kernel_check_childs(u, pattern_vid);
    }

    void update(graphchi_vertex<VertexDataType, EdgeDataType> &vertex, graphchi_context &gcontext) {
        
        if (scheduler) gcontext.scheduler->remove_tasks(vertex.id(), vertex.id());
        
        if (gcontext.iteration == 0) {
            vertex.set_data(rand() % 5);
            if (scheduler)  gcontext.scheduler->add_task(rand() % 5);
        }
        
        if (gcontext.iteration == 1) {
            if (match_set_->indicator_->get_bit(vertex.vertexid)) ;
            else {
                for(size_t pattern_i = 0; pattern_i < pattern_->num_vertexes_; pattern_i++) {
                    if(vertex.get_data() == pattern_->vertexes_[pattern_i]->label_) {
                        match_set_->indicator_->set_bit(vertex.vertexid);
                        match_set_->sim_sets_[vertex.vertexid]->set_bit(pattern_i);
                        //printf("%d match %d, label: %d\n ", vertex.vertexid, pattern_i, pattern_->vertexes_[pattern_i]->label_);
                    }
                }
            }
        }
        /* Check childs */
        if (gcontext.iteration > 1 ) {
            for(size_t i = 0; i< pattern_->num_vertexes_; i++)
            {
                auto pattern_vid = pattern_->vertexes_[i]->vid_;
                if(match_set_->sim_sets_[i]->get_bit(pattern_vid))
                {
                    if(!CheckU(vertex, (unsigned)pattern_vid))
                    {
                        match_set_->sim_sets_[i]->rm_bit(pattern_vid);
                        for(int j=0; j < vertex.num_edges(); j++) {
                        if (scheduler) gcontext.scheduler->add_task(vertex.edge(j)->vertex_id(), true);
                        converged = false;

                        }
                    }
                }

            }
        }


    }    
    /**
     * Called before an iteration starts.
     */
    void before_iteration(int iteration, graphchi_context &info) {
        iterationcount++;
        converged = iteration > 0;
    }
    
    /**
     * Called after an iteration has finished.
     */
    void after_iteration(int iteration, graphchi_context &ginfo) {
        if (converged) {
            std::cout << "Converged!" << std::endl;
            ginfo.set_last_iteration(iteration);
        }
    }
    
    /**
     * Called before an execution interval is started.
     */
    void before_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
    }
    
    /**
     * Called after an execution interval has finished.
     */
    void after_exec_interval(vid_t window_st, vid_t window_en, graphchi_context &ginfo) {        
    }
    
};

int main(int argc, const char ** argv) {
    /* GraphChi initialization will read the command line 
     arguments and the configuration file. */
    graphchi_init(argc, argv);

    /* Metrics object for keeping track of performance counters
     and other information. Currently required. */
    metrics m("connected-components");
    
    /* Basic arguments for application */
    std::string filename = get_option_string("file");  // Base filename
    int niters           = get_option_int("niters", 1000); // Number of iterations (max)
    scheduler            = get_option_int("scheduler", false);
    
    /* Process input file - if not already preprocessed */
    int nshards             = (int) convert_if_notexists<EdgeDataType>(filename, get_option_string("nshards", "auto"));
    
    int n = 5;
    printf("InitPatternNClique\n");
    ImmutableCSR * pattern = InitPatternNClique(n);
    MatchSet * matchset = new MatchSet(2147483647, n, true);
    printf("ShowImmutableCSR");
    ShowImmutableCSR(*pattern);
    
    printf("#################################\n");
    
    if (get_option_int("onlyresult", 0) == 0) {
        /* Run */
        SimPragram program;
        program.match_set_ = matchset;
        program.pattern_ = pattern;

        graphchi_engine<VertexDataType, EdgeDataType> engine(filename, nshards, scheduler, m); 
        engine.run(program, niters);
    }
    
    /* Run analysis of the connected components  (output is written to a file) */
    m.start_time("label-analysis");
    
    analyze_labels<vid_t>(filename);
    m.stop_time("label-analysis");
    
    /* Report execution metrics */
    metrics_report(m);
    return 0;
}

