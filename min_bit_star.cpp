// minimum viable bit star
// thomas ck fuller.  
// "an expert is someone whos made all the mistakes that are possible"
#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <cmath>

// i think i want to do a header file that contains my reading in code or something
// there is probably a smart way to keep everything organized

class BITSTAR{
    public:
    private:
        // putting debugging utilities in the "private" part of the class structure
    public:
        // major params for BIT* 

        // Boundaries
        float xMax = 10.0f; 
        float yMax = 10.0f; 

        // Sampling Parameters
        int m = 50;             // number of samples

        // RGG Parameters
        int kBitStar = 10; 
        float rBitStar = 3.0f;

        // Calculate L2 Norm
        float calculate_L2(float x1, float y1, float x2, float y2){
            float L_2 = sqrt( pow((x2-x1),2) + pow((y2-y1),2) ); 
            return L_2;
        };
        // all of our data structures for the BIT* implementation go here.
        struct state { // i guess i just add more members there if i need it
            float x = INFINITY;                 // valgrind didnt like when i didn't initialize stuff
            float y = INFINITY;
            float f = INFINITY;
            float gT = INFINITY;
            float hHat = INFINITY;
            float gHat = INFINITY;              // structure equality just checks a subset of these
            float r = INFINITY;                // distance for nearest neighbor search 
            // associated heuristic information goes here?
        };
        struct edge {
            state source_state;
            state target_state; 
            float f = INFINITY;
            float cHat;
        };

        float fCalculateDist(state v, state goal){
            float xv, yv, xg, yg;                       // idk why but this feels smart? but also dumb
            xv = v.x; yv = v.y;
            xg = goal.x; yg = goal.y;
            float dist = calculate_L2(xv,yv,xg,yg);
            return dist;
        };
        //struct tree {                 
        //    std::vector<state> V;
        //    std::vector<edge> E;
        //};

        // priority functions for the vertex queue
        struct neighborQueueSort {
            bool operator() (state const& s1, state const& s2){
                return s1.r > s2.r;                                             // this ordering produces a min queue, which is what I need
            }
        };
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

        // i love typedef
        typedef std::priority_queue<state, std::vector<state>, vertexQueueSort> vertexQueueType;             // 3.0
        typedef std::priority_queue<state, std::vector<state>, neighborQueueSort> neighborQueueType;             // 3.0
        typedef std::priority_queue<edge, std::vector<edge>, edgeQueueSort> edgeQueueType;                 // 3.1
        typedef std::vector<edge> edgeVector;
        typedef std::vector<state> stateVector;

        // something from stackoverflow that i hope clears my queue
        template<class Q>           // dont really understand how template works but i should probably learn it
        void clearQueue(Q& q) {
            q = Q();
        }

        // BIT_STAR
        void BIT_STAR(state start, state goal){
            // preprocessing the start vertex
            start.gT = 0.0f; 
            start.hHat = fCalculateDist(start,goal);
            start.f = start.gT + start.hHat;                            // do i need this every time?`
            //std::cout << "start.f: " << start.f << std::endl;

            // preprocessing for the goal state
            goal.hHat = 0.0f;

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
            Vunexpnd.push_back(start);                                                       // 4.1 
            std::vector<state> Xnew;                                                        // 4.2
            Xnew.push_back(goal);                                                           // 4.2 this maybe only working because the goal right now is a single point.

            float ci = INFINITY;                                                            // 5.0  This feels like cheating, since on the first pass it will definitely be equal to infinity

            // stuff for lines 7-12
            stateVector Xsampling;
            stateVector Xreuse;
            
            // REPEAT                                                                       // 6.0
            while (true) {
                std::cout << "Qe.size() == " << Qe.size() << " Qv.size() == " << Qv.size() << std::endl;
                if ((Qe.size() == 0) && (Qv.size() == 0)){                                      // 7.0
                    std::cout << "The part where I sample stuff " << std::endl;
                    // we are going to write the pruning function last i think
                    // Xreuse = Prune()                                                         // 8.0 (not implemented yet)
                    Xsampling = Sample(m, start, goal, ci) ;                                    // 9.0; 
                    Xnew = Xsampling;                                                           // 10.0 (need to union with reuse states when the time comes)
                    // i kind of like micromanaging data structures but also i hate it
                     



                    std::exit(0); // causes 120 bytes to be lost on exit
                    // 12.0 
                };
                while  ((Qv.size() > 0) &  (fV_BestQueueValue(Qv) <= fE_BestQueueValue(Qe)) ){   // 13.0
                //  on the first iteration, i shoudl be in here, since the vertex queue is empty
                    // EXPAND NEXT VERTEX
                    ExpandNextVertex(Qv, Qe, ci, Vunexpnd, Xunconn, start, goal,V,E);                                               //  14.0
                    std::cout << " || || | ||" << std::endl;
                    std::cout << "\tBEST VERTEX VALUE " << fV_BestQueueValue(Qv);
                    std::cout << "\tBEST EDGE VALUE " << fE_BestQueueValue(Qe) << std::endl;
                    std::cout << "\tNUMBER OF VERTICES IN V " << V.size()<< std::endl;
                    std::cout << "\tNUMBER OF EDGES IN E " << E.size() << std::endl;

                    std::cout << "\tVertex Queue Size " << Qv.size() << std::endl;
                    std::cout << "\tEDGE Queue Size " << Qe.size() << std::endl;
                    std::cout << "\tVunexpnd Size " << Vunexpnd.size() << std::endl;  
                    std::cout << " || || | ||" << std::endl;
                };
                std::cout << " I think its failing when it tries to pop into an empty edge queue "  << std::endl;

                if (!Qe.empty()) {    // it only makes sense to pop from a non empty queue
                    edge currentEdge = E_PopBestInQueue(Qe); 
                    std::cout << "Trying to understand whats in the edge queue" << std::endl;
                    std::cout << "source state: "; 
                    std::cout << currentEdge.source_state.x << " " << currentEdge.source_state.y << std::endl;
                    std::cout << "target state: ";
                    std::cout << currentEdge.target_state.x << " " << currentEdge.target_state.y << std::endl;
                } else {
                    // how to empty a queue                    
                    clearQueue(Qe); clearQueue(Qv);                                                                     // 34.0
                    std::cout << "is the vertex queue empty? " << Qv.empty() << std::endl;
                    std::cout << "is the edge queue empty? " << Qe.empty() << std::endl;

                } 
                
                //std::cout << "diagnostic message to let you know i am stuck in this while loop" <<std::endl;
                //std::cout << " how big is my edge queue " << Qe.size() << std::endl;
            };//UNTIL STOP;
            // RETURN T;
        }; // MAIN BIT_STAR END

        // EXPAND NEXT VERTEX (Algorithm 2)
        void ExpandNextVertex(vertexQueueType& Qv, edgeQueueType& Qe, 
                              float& ci, stateVector& Vunexpnd, stateVector& Xunconn,
                              state& start, state& goal, stateVector& V, edgeVector& E){                                       // there has to be a cleaner way
            std::cout << "Expanding next vertex. Cost = " << ci << std::endl;
            stateVector Xnear, XnearAndUnconnected, Vnear; 
            state Vmin = sV_PopBestInQueue(Qv);                                                 // A2.1
            float gHat; 
            float cHat; 
            float hHat;
            //if Vmin in Vunexpanded 
            bool vminInVunexpnd = b_VertexIsIn(Vmin, Vunexpnd);

            // check the samples
            if (vminInVunexpnd) {                                                               // A2.2
                std::cout << "vertex is IN!" << std::endl;
                Xnear = Near(Xunconn, Vmin, rBitStar, kBitStar);                                // A2.3
                std::cout << "X near: " << Xnear.size(); 
            }
            else {                                                                              // A2.4
                std::cout << "vertex is OUT!" << std::endl;
                // scary since i wont be able to test this right away
                XnearAndUnconnected = SetIntersection(Xnear, Xunconn);                          // A2.5
                Xnear = Near(XnearAndUnconnected, Vmin, rBitStar, kBitStar);                    // A2.5
            }
            for (auto &i : Xnear) { 
                gHat = fCalculateDist(i, start); i.gHat = gHat;                                 // A2.6
                cHat = fCalculateDist(i, Vmin);                                                 // A2.6 
                hHat = fCalculateDist(i, goal); i.hHat = hHat;                                  // A2.6
                if (gHat + cHat + hHat < ci) {                                                  // A2.6
                    edge EdgeToPush;                                                            // A2.6
                    EdgeToPush.source_state = Vmin; EdgeToPush.target_state = i;                // A2.6
                    //EdgeToPush.f = gHat + cHat + hHat;
                    EdgeToPush.cHat = cHat;
                    Qe.push(EdgeToPush);                                                        // A2.6
                    std::cout << "enqueing edge based on Xnear" << std::endl;
                };
            } // A2.6  for loop

            // check the motion tree
            if (vminInVunexpnd) {                                                               // A2.7
                std::cout << "in the second part of the loop" << std::endl;
                Vnear = Near(V, Vmin, rBitStar, kBitStar);                                      // A2.8
                std::cout << Vnear.size() << std::endl;
                for (auto &i : Vnear) {
                    gHat = fCalculateDist(i, start); i.gHat = gHat;                             // stuff to make 2.9 possible
                    cHat = fCalculateDist(i, Vmin);                                             // 
                    hHat = fCalculateDist(i, goal); i.hHat = hHat;                              // 
                    edge EdgeToPush;                                                            // 
                    EdgeToPush.source_state = Vmin; EdgeToPush.target_state = i;                // 
                    if (gHat + cHat + hHat < ci && !b_EdgeIsIn(EdgeToPush, E)) {                // A2.9
                        EdgeToPush.f = gHat + cHat + hHat;
                        EdgeToPush.cHat = cHat; 
                        Qe.push(EdgeToPush);
                        std::cout << "enqueing edge based on Vnear" << std::endl;
                    }; 
                }// A2.9 for loop
            };
            // remove Vmin from the "unexpanded list"                                           
            removeStateFromSet(Vmin, Vunexpnd);                                                 // A2.10
        }; 

        stateVector SetIntersection(stateVector set1, stateVector set2) {
            // inputs: takes two state vectors and outputs the common states
            // this is O(n x m), but idk how to do it faster
            // yi suggested a hash table, but I'm just going to use this for now.
            stateVector commonStates; 
            for (auto &i : set1) {
                    if (b_VertexIsIn(i, set2)) commonStates.push_back(i); 
            }; 
            return commonStates;
        };
        stateVector Sample(int m, state start, state goal,  float ci){
            // inputs m: int: how many samples
            // state start: start, state goal: goal state
            // float ci: the cost of the informative region, calculated by l2(goal, rand) + l2(start,rand)
            // can softlock if the cost ci is too low
            stateVector sampledStates;
            // sample states and if they are within 
            float xRand;
            float yRand;
            float xs, ys;       // this is so much effor to test other types of systems!!
            float xg, yg;
            float tmpG, tmpH;
            int counter = 0;
            xs = start.x; ys = start.y;         // am i rushing?
            xg = goal.x; yg = goal.y;
            std::cout << " X|Y|C" << std::endl;
            while (counter < m) {

                xRand = (float(rand()) / float(RAND_MAX))*xMax;     // trims the samples to the boundary
                yRand = (float(rand()) / float(RAND_MAX))*yMax;
                tmpG = calculate_L2(xs,ys,xRand,yRand);             // calculating the ellispe
                tmpH = calculate_L2(xg,yg,xRand,yRand); 
                if ((tmpG + tmpH) < ci) {
                    // theres going to be collision checking based on the input MAP that needs to go here
                    // the FLOOR function that i studied was written in notebook 2 page 105
                    state stateToAdd;
                    stateToAdd.x = xRand; stateToAdd.y = yRand;
                    sampledStates.push_back(stateToAdd);
                    counter++;      // only increment the counter if the sample was valid
                
                    std::cout << stateToAdd.x << " | " << stateToAdd.y << " | " << tmpG+tmpH <<  std::endl;
                }
            }; 

            return sampledStates;
        }; 
        stateVector Near(stateVector StatesToCheck, state vertexToCheck,  float rggRadius, int numberOfNeighbors) {
            // I need to pick the RGG radius or the KNN
            // inputs:  statevector: to be searched
            //          state: measuring nearness from this point
            //          float: a radius to say yes this is within my rgg radius
            //          int: the number of nearest states 
            //  outputs: the list of numberOfNeighbors nearest states
            // i'm just praying that the subset of "states to check" doesn't bog down the alg.
            stateVector NearestNeighbors;
            neighborQueueType NeighborQueue;                                // sort on insertion     
            int neighborCounter = 0;
            for (auto &i : StatesToCheck) {
                i.r = fCalculateDist(i, vertexToCheck);                 // i = thing i'm iterating
                NeighborQueue.push(i);
            };

            // when things get put on the neighbor queue they are already sorted, so 
            // i THINK (this may not be happening), that its safe to increment the neighbor counter every time
            while (!NeighborQueue.empty() && neighborCounter < numberOfNeighbors) {
                state TopState = NeighborQueue.top();
                NeighborQueue.pop();
                if (TopState.r < rggRadius && !b_VerticesAreEqual(TopState, vertexToCheck)) {
                    NearestNeighbors.push_back(TopState);
                    neighborCounter++;
                }
            };
            return NearestNeighbors;
        }; // Near

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
            std::cout << "pop best in vertex queue " << std::endl; 
            state Top = Qv.top();
            std::cout << Top.x << std::endl;;
            Qv.pop();
            return Top;
        };
        edge E_PopBestInQueue(edgeQueueType& Qe) {
            // there has to be a smart way to do this
            std::cout << "pop best in edge queue" << std::endl;
            edge Top = Qe.top();
            Qe.pop();
            return Top;
        };

        bool b_VertexIsIn(state V, stateVector vectorSet){ // i have it with reference, but since its retunring i think its safe to not do that
            // input: a state vertex and an associate set
            // output: return true if the state vector contains that thing, else false
            bool verticesAreEqual;
            for (auto &r : vectorSet) {
                //std::cout << "trying this" << std::endl;
                verticesAreEqual = b_VerticesAreEqual(r, V);
                if (verticesAreEqual) {
                    //std::cout << "set membership is shown" << std::endl;
                    return verticesAreEqual;
                };
            };
            return false; 
        };

        bool b_VerticesAreEqual(const state& V1, const state& V2) {
            // input is two states: output is whether or not every member of their structure is equal
            bool verticesAreEqual = ((V1.x == V2.x) && 
                                    (V1.y == V2.y) &&
                                    (V1.f == V2.f) &&
                                    (V1.gT == V2.gT) &&
                                    (V1.hHat == V2.hHat));
            return verticesAreEqual; 
        };
        
        bool b_EdgesAreEqual(edge E1, edge E2) {
            bool edgesAreEqual = b_VerticesAreEqual(E1.source_state, E2.source_state) &&
                                 b_VerticesAreEqual(E1.target_state, E2.target_state);
            return edgesAreEqual;
        };
        bool b_EdgeIsIn (edge E, edgeVector edgeSet) {
            // input is an edge and an edge vector
            // output: whether or not the edge is a member of the corresponding edge set
            // equality of the edge is measured if both states are the same 
            bool edgeIsIn; 
            for (auto &r : edgeSet) {
                edgeIsIn = b_EdgesAreEqual(r, E);
                if (edgeIsIn) {
                    return edgeIsIn; 
                }
            };
            return false;
        };

        bool removeStateFromSet(state Vertex, stateVector& VectorList) {
            // inputs: a Vertex and a stateVector
            // effect: it removes that state from the stateVector VectorList
            // if we dont see that effect globally, we can trace the references back
            // the VectorList is passed by reference to modify the ACTUAL list that you are working with.
            int stateIdx = 0;
            for (auto &i : VectorList) {
                // asdfadsf
                std::cout << "in removeStateFromSet" << std::endl;
                if (b_VerticesAreEqual(i, Vertex)){
                    VectorList.erase(VectorList.begin() + stateIdx);
                    return true;    // so we are hoping that multiple copies dont get in here
                }
                stateIdx++; 
            };
            return false;
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

    float calculateGT(state& stateOfInterest, state& start, stateVector& V, edgeVector& E) {
        // inputs: state: a state you want to find the cost to get to it given the current tree
        //         edgeVector: a list of all the edges in the current motion tree
        //         stateVector: a list of the states in the current motion tree
        // effects: this function attempts to search the motion tree until the root is reached, 
        // page 111 of notebook 2 shows the tree that this wokrs on.
        float gT = 0.0f; 
        bool stateFound = false;
        // this is hard for me to visualize, how to instruct the computer what to do. 
        state tmpSourceState;
        while (!b_VerticesAreEqual(stateOfInterest, start)){         // hopefully stops my tree traversal when i hit the start state
            std::cout << "im in the while loops stuck for good" << std::endl;
            stateFound = false;                                // if no state is found while searching the edges, we need to return gT = INFINITY
            for (auto &i : E) {
                if(b_VerticesAreEqual(stateOfInterest, i.target_state)){
                    stateFound = true;        
                    gT +=i.cHat;                                    // i think this works since states in the motion tree have non infinite costs
                    std::cout << "gT value: " << gT << std::endl;  // it never looks 
                    stateOfInterest = i.source_state;
                }; 
            };
            if (!stateFound) return INFINITY;                       // need this for if hte tree isnt connected
        }
       return gT;
    };

    stateVector Append(stateVector V1, stateVector V2){
        stateVector tmp1, tmp2;
        if (V1.size() > V2.size() ) {
            tmp1 = V1; tmp2 = V2;
        }
        else {
            tmp1 = V2; tmp2 = V1; 
        };
        tmp1.insert(tmp1.end(), tmp2.begin(), tmp2.end());
        return tmp1;
    };

    // UNIT TESTS OF EACH PART 
    bool bTest_b_VerticesAreEqual() {
        state V1, V2;  
        V1.x = 1.0f; V2.x = 1.0f; 
        V1.y = 0.0f; V2.y = 0.0f; 
        V1.f = 3.14f; V2.f = 3.14;  
        V1.gT = INFINITY; V2.gT = INFINITY; 
        V1.hHat = INFINITY; V2.hHat = INFINITY;
        bool areEqual = b_VerticesAreEqual(V1, V2);
        return areEqual;
    };
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

    bool bTest_Near() {
        stateVector NeighborsCheck; 
        stateVector StatesToCheck; 
        state vertexToCheck;
        float rggRadius = 3.0f;
        int numberOfNeighbors = 2; 
        // i'll test 3 states, and see if just two appear in the neighbor list
        state s0, s1, s2, s3, s4, s5, s6;
        s0.x = 0.0f; s0.y = 0.0f;
        s1.x = 1.0; s1.y = 1.0; 
        s2.x = -1.5; s2.y = -1.5; 
        s3.x = 3.0; s3.y = 3.0;
        s4.x = -1.6; s4.y = 1.6;
        s5.x = -1.0; s5.y = 1.0;
        s6.x = 0.0f; s6.y = 0.0f;
        StatesToCheck.push_back(s5); StatesToCheck.push_back(s6);   
        StatesToCheck.push_back(s1);StatesToCheck.push_back(s4);  StatesToCheck.push_back(s3); StatesToCheck.push_back(s2);
        NeighborsCheck = Near(StatesToCheck, s0, rggRadius, numberOfNeighbors); 
        std::cout << NeighborsCheck.size() << std::endl; 
        std::cout << "X | Y | r " << std::endl;

        for (auto &i : NeighborsCheck) {
            std::cout << i.x << " | " << i.y << " | " << i.r << std::endl;
        };
        if (NeighborsCheck.size() == 2) return true;            // not a perfect check but better than nothing
        

        return false;
    };
    bool bTest_SetIntersection(){
        // if s8 and s1 are apparently added to commonStates, we succeed
        stateVector commonStates;
        stateVector set1, set2;   
        state s1, s2, s3, s4, s5, s6, s7, s8; 
        s1.x = 1.0f; s1.y = 1.0f;
        s2.x = 2.0f; s2.y = 2.0f;
        s3.x = 3.0f; s3.y = 3.0f;
        s4.x = 4.0f; s4.y = 4.0f;
        s5.x = 5.0f; s5.y = 5.0f;
        s6.x = 6.0f; s6.y = 6.0f;
        s7.x = 7.0f; s7.y = 7.0f;
        s8.x = 8.0f; s8.y = 8.0f;
        set1.push_back(s8); set2.push_back(s8);
        set1.push_back(s1); set2.push_back(s2); 
        set1.push_back(s3); set2.push_back(s4); 
        set1.push_back(s5); set2.push_back(s6); 
        set1.push_back(s7); set2.push_back(s1); 
        commonStates = SetIntersection(set1, set2);
        if (commonStates.size() == 2) return true;
        return false;
    };
    
    bool bTest_b_EdgesAreEqual() {
        // tests edge equality and edge membership
        // this also tests edge membership
        edgeVector edgeList; 
        edge E1, E2, E3;
        state E1s, E1t, E2s, E2t, E3s, E3t; 
        E1s.x = 0.0f; E1s.y = 0.0f;
        E1t.x = 0.0f; E1t.y = 0.0f;

        E2s.x = 0.0f; E2s.y = 0.0f;
        E2t.x = 0.0f; E2t.y = 0.0f;
        
        E3s.x = 0.0f; E3s.y = 0.0f;
        E3t.x = 1.0f; E3t.y = 1.0f;

        E1.source_state = E1s; E1.target_state = E1t;
        E2.source_state = E2s; E2.target_state = E2t;
    
        edgeList.push_back(E1); edgeList.push_back(E2); edgeList.push_back(E3); 
        bool isTheEdgeIn = b_EdgeIsIn(E2, edgeList); 
        
        if (b_EdgesAreEqual(E1, E2) && !b_EdgesAreEqual(E1,E3) && isTheEdgeIn) return true;
        return false;


    };

    bool bTest_removeStateFromSet(){
        // testing state removal from vector, kinda wish i messed around with hash table now, 
        // but i'm committed
        state s1, s2, s3, s4; 
        stateVector stateList;
        s1.x = 1.0f; s1.y = 1.0;
        s2.x = 1.0f; s2.y = 2.0;
        s3.x = 3.0f; s3.y = 3.0;
        s4.x = 4.0f; s4.y = 4.0;
        stateList.push_back(s1); stateList.push_back(s2); 
        stateList.push_back(s3); stateList.push_back(s4); 
        removeStateFromSet(s4, stateList); 
        print_state_vector(stateList, "stateList");
        // read the printout to make sure
        if (stateList.size() == 3 && stateList[2].x == 3.0) return true;
        return false;
    };

    bool bTest_calculateGT() {
        // i kind of need gT to be calculated for a functioning algorithm
        edge E0, E1, E2, E3, E4, E5;
        state start, s1, s2, s3, s4, s5;
        stateVector V; edgeVector E;
        start.x = 0.0f; start.y = 0.0f;
        s1.x = 1.0f; s1.y = 1.0;
        s2.x = 1.0f; s2.y = 2.0;
        s3.x = 3.0f; s3.y = 3.0;
        s4.x = 4.0f; s4.y = 4.0;
        s5.x = 5.0f; s5.y = 5.0; 
        E0.source_state = start; E0.target_state = s1;   E0.cHat = 4.0f;
        E1.source_state = s1; E1.target_state = s2;      E1.cHat = 3.0f; 
        E2.source_state = s1; E2.target_state = s3;      E2.cHat = 2.0f;
        E3.source_state = s3; E3.target_state = s4;      E3.cHat = 5.0f;
        E4.source_state = s3; E4.target_state = s5;      E4.cHat = 7.0f; 
        V.push_back(start);
        V.push_back(s1); V.push_back(s2); V.push_back(s3); V.push_back(s4); V.push_back(s5);
        E.push_back(E0); 
        E.push_back(E1); E.push_back(E2); E.push_back(E3); E.push_back(E4);
        float costGivenCurrentTree = calculateGT(s5, start, V, E);
        std::cout<< "Cost Given Current test tree: " << costGivenCurrentTree << std::endl; 
        if (costGivenCurrentTree == 13.0f) return true;
        return false;

    };
    
    bool bTest_Sample() {
        // not sure how I'm going to be able to do this one

        state start; state goal;
        start.x = 1.0; start.y = 1.0;
        goal.x = 8.0; goal.y = 8.0;
        int m = 50;
        float ci = 12.0f;
        Sample(50, start, goal, ci);
        return true; 
    };

    bool bTest_Append() {
        stateVector Vec1, Vec2; 
        stateVector resultVector;
        state s1, s2, s3, s4, s5;
        s1.x = 1.0f; s1.y = 1.0;
        s2.x = 1.0f; s2.y = 2.0;
        s3.x = 3.0f; s3.y = 3.0;
        s4.x = 4.0f; s4.y = 4.0;
        s5.x = 5.0f; s5.y = 5.0; 
        Vec1.push_back(s1); Vec1.push_back(s2);
        Vec2.push_back(s3); Vec2.push_back(s4);
        Vec2.push_back(s5); 
        
        resultVector = Append(Vec1, Vec2);
        if (resultVector.size() == 5) return true;
        return false;
 
    };
    bool unit_test() {
        bool bqvv_works = bTest_fV_BestQueueValue();
        bool bqve_works = bTest_fE_BestQueueValue();
        bool bvae_works = bTest_b_VerticesAreEqual();
        bool near_works = bTest_Near();
        bool sint_works = bTest_SetIntersection();
        bool edge_equal_works = bTest_b_EdgesAreEqual();
        bool state_V_removal_works = bTest_removeStateFromSet();
        bool tree_traversal_works = bTest_calculateGT();
        bool sampling_works = bTest_Sample();
        bool appending_works = bTest_Append();

        std::cout << std::endl;
        std::cout << "UNIT TEST RESULTS: " << std::endl;
    
        std::cout << "\tBest Vertex Queue Value: " << bqvv_works << std::endl;  
        std::cout << "\tBest Edge Queue Value: " << bqve_works << std::endl;  
        std::cout << "\tVertex Equality: " << bvae_works << std::endl;  
        std::cout << "\tNear() Function: " << near_works << std::endl;  
        std::cout << "\tSet Intersection / Vertex Membership: " << sint_works << std::endl;
        std::cout << "\tEdge Equality: " << edge_equal_works << std::endl;
        std::cout << "\tState Removal From V Vector: " << state_V_removal_works << std::endl;
        std::cout << "\tTree Traversal Works: " << tree_traversal_works << std::endl;
        std::cout << "\tAppending Works: " << appending_works << std::endl;
        
        return bqvv_works && bqve_works && 
               bvae_works && near_works && sint_works && 
               edge_equal_works && state_V_removal_works &&
               tree_traversal_works && appending_works;
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
        
    std::cout << std::endl;
    std::cout << "\tTEST WORKED IF 1 == " << test_works << std::endl;
    std::cout << std::endl;
    if (test_works) {
        BS.BIT_STAR({2.1, 3.0f}, {9.0f, 9.0f}); // start (x,y) = (1.1, 3.0), goal (9,9)
    };

    // unit testing my stuff
    return 0;
};
