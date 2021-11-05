// minimum viable bit star
// thomas ck fuller.  
// "an expert is someone whos made all the mistakes that are possible"
#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <cmath>
#include <limits>

class BITSTAR{
    public:
    private:
        // putting debugging utilities in the "private" part of the class structure
    public:
        // all our math helper functions going to go at the top
        float calculate_L2(float x1, float y1, float x2, float y2){
            float L_2 = sqrt( pow((x2-x1),2) + pow((y2-y1),2) ); 
            return L_2;
        };
        // all of our data structures for the BIT* implementation go here.
        struct state { // i guess i just add more members there if i need it
            float x;
            float y;
            float f;
            float gT = INFINITY;
            float hhat = INFINITY;
            // associated heuristic information goes here?
        };
        struct edge {
            state source_state;
            state target_state; 
            float f = 0.13131;
            float chat;
        };

        float fCalculateHhat(state v, state goal){
            float xv, yv, xg, yg;                       // idk why but this feels smart? but also dumb
            xv = v.x; yv = v.y;
            xg = goal.x; yg = goal.y;
            float hhat = calculate_L2(xv,yv,xg,yg);
            return hhat;
        };
        //struct tree {                 
        //    std::vector<state> V;
        //    std::vector<edge> E;
        //};

        // priority functions for the vertex queue
        struct vertexQueueSort {
            bool operator() (state const& s1, state const& s2){
                return s1.f > s2.f;                                             // this ordering produces a min queue, which is what I need
            }
        };
        struct edgeQueueSort {
            bool operator() (edge const& e1, edge const& e2){                   // sort occurs on an edge-based H E U R I S T I C 
                return e1.f > e2.f;
            }
        };

        // i fuck fuck fucking love typedef
        typedef std::priority_queue<state, std::vector<state>, vertexQueueSort> vertexQueueType;             // 3.0
        typedef std::priority_queue<edge, std::vector<edge>, edgeQueueSort> edgeQueueType;                 // 3.1
        typedef std::vector<edge> edgeVector;
        typedef std::vector<state> stateVector;

        // BIT_STAR
        void BIT_STAR(state start, state goal){
            // preprocessing the start vertex
            start.gT = 0.0f; 
            start.hhat = fCalculateHhat(start,goal);
            start.f = start.gT + start.hhat;                            // do i need this every time?`
            std::cout << "start.f: " << start.f << std::endl;

            // preprocessing for the goal state
            goal.hhat = 0.0f;

            // algorithm begins as in the paper
            stateVector V;                                                                  // 1.0 
            V.push_back(start);                                                             // 1.0
            edgeVector E;                                                                   // 1.2

            stateVector Xunconn;                                                            // 2.0
            Xunconn.push_back(goal);                                                        // 2.0

            vertexQueueType Qv;                                                             // 3.0
            Qv.push(start);                                                                 // 3.0
            edgeQueueType Qe;                                                               // 3.1

            stateVector Vsoln;                                                              // 4.0
            Vsoln.push_back(goal);                                                          // 4.0 
            std::vector<state> Vunexpnd;                                                    // 4.1
            Vunexpnd.push_back(goal);                                                       // 4.1 
            std::vector<state> Xnew;                                                        // 4.2
            Xnew.push_back(goal);                                                           // 4.2 this maybe only working because the goal right now is a single point.

            float ci = INFINITY;                                                            // 5.0  This feels like cheating, since on the first pass it will definitely be equal to infinity

            // REPEAT                                                                       // 6.0
            while (true) {
                if ((Qe.size() == 0) && (Qv.size() == 0)){                                      // 7.0
                    std::cout << "testing if I got into the lines 7-12" << std::endl;
                    std::cout << "lines 7 thru 12 don't happen the first iteration" << std::endl;
                    std::exit(0); // causes 120 bytes to be lost on exit
                    // 12.0 
                };
                while  ((Qv.size() > 0) &  (fV_BestQueueValue(Qv) <= fE_BestQueueValue(Qe)) ){   // 13.0
                //  on the first iteration, i shoudl be in here, since the vertex queue is empty
                    // EXPAND NEXT VERTEX
                    ExpandNextVertex(Qv, Qe, ci, Vunexpnd);                                               //  14.0
                    std::cout << "BEST VERTEX VALUE " << fV_BestQueueValue(Qv);
                    std::cout << " BEST EDGE VALUE " << fE_BestQueueValue(Qe) << std::endl;
                    std::cout << "  Vertex Queue Size " << Qv.size();
                    std::cout << "  EDGE Queue Size " << Qe.size();
                };
            };//UNTIL STOP;
            // RETURN T;
        }; // MAIN BIT_STAR END
        void ExpandNextVertex(vertexQueueType& Qv, edgeQueueType& Qe, float& ci, stateVector& Vunexpnd){
            std::cout << "were expanding next vertex " << ci << std::endl;
            state Vmin = sV_PopBestInQueue(Qv);                                                 // A2.1
            //if Vmin in Vunexpanded 
            if (b_VertexIsIn(Vmin, Vunexpnd)) std::cout << "vertex is IN!" << std::endl;
            std::cout << "size of the queue in ExpandNextVertex" << Qv.size() << std::endl;
        }; 
        float fV_BestQueueValue(vertexQueueType& Qv) {              
            if (Qv.size() == 0) return INFINITY;
            return Qv.top().f;
        };
        float fE_BestQueueValue(edgeQueueType& Qe){                
            if (Qe.size() == 0) return INFINITY;
            return Qe.top().f;
        };
        state sV_PopBestInQueue(vertexQueueType& Qv){
            // will i regret not writing a unit test for this?
            std::cout << "pop best in queue " << std::endl; 
            state Top = Qv.top();
            std::cout << Top.x << std::endl;;
            Qv.pop();
            return Top;
        };
        bool b_VertexIsIn(state& V, stateVector& vectorSet){
            // input: a state vertex and an associate set
            // output: return true if the state vector contains that thing, else false
            for (auto &r : vectorSet) {
                std::cout << "trying this" << std::endl;
            };
            return true; 
        };
        bool b_VerticesAreEqual(const state& V1, const state& V2) {
            // input is two states: output is whether or not every member of their structure is equal
            return true; 
        };
        // I should write simple test functions as I go along
        void print_state_vector(std::vector<state> state_vector, std::string vector_name) {
            int i=0; 
            std::cout << vector_name << std::endl;
            for (auto &r : state_vector){
                std::cout << "iteration: " << i  << " member x, " << r.x << " member y " << r.y << std::endl; 
                i++;
            };
        };

    // UNIT TESTS OF EACH PART 
    bool bTest_fV_BestQueueValue() {              
        float upperValue = 10.0f; 
        float middleValue = 9.0f;
        float lowerValue = 8.0f;
        bool check_inf, check_val; 
        std::priority_queue<state, std::vector<state>, vertexQueueSort> Qv;            

        if (fV_BestQueueValue(Qv) == INFINITY) check_inf = true;
        state v1 = {0.0f, 0.0f, upperValue};
        Qv.push(v1);

        state v2 = {0.0f, 0.0f, lowerValue};
        Qv.push(v2);

        state v3 = {0.0f, 0.0f, middleValue};
        Qv.push(v3);

        float a = fV_BestQueueValue(Qv);
        if (fV_BestQueueValue(Qv) == lowerValue) check_val = true;

        if (check_val & check_inf) return true;
        return false;
    };

    bool bTest_fE_BestQueueValue() {
        float upperValue = 420.0f;
        float middleValue = 100.0f;
        float lowerValue = 69.0f;
        bool check_inf, check_val;

        std::priority_queue<edge, std::vector<edge>, edgeQueueSort> Qe;                 // 3.1
        if (fE_BestQueueValue(Qe) == INFINITY) check_inf = true;
        edge e1; edge e2; edge e3;
        e1.f = upperValue; 
        e2.f = lowerValue;
        e3.f = middleValue;
        Qe.push(e1); Qe.push(e2); Qe.push(e3); 
        float a = fE_BestQueueValue(Qe);
        if (fE_BestQueueValue(Qe) == lowerValue) {
            check_val = true;
        }
        if (check_val & check_inf) return true;
        return false;
    };

    bool unit_test() {
        bool bqvv_works = bTest_fV_BestQueueValue();
        bool bqve_works = bTest_fE_BestQueueValue();
        return bqvv_works && bqve_works;
    };

//// BIT STARS ENDING BRACE DO NOT TOUCH
};///DONT TOUCH
//// DONTTTT TOUCH

int main () {
    // should i invest the time into writing a test harness? the answer is almost definitely yes

    // the main program
    std::cout << "minimum viable bit star" << std::endl;
    BITSTAR BS;
    bool test_works;
    test_works = BS.unit_test(); 
        
    std::cout << "did the test work ? : " << test_works << std::endl;
    if (test_works) {
        BS.BIT_STAR({1.1, 3.0f}, {9.0f, 9.0f});
    };

    // unit testing my stuff
    return 0;
};
