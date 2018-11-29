#include <bits/stdc++.h>
#include <time.h>
using namespace std;
// Format checker just assumes you have Alarm.bif and Solved_Alarm.bif (your file) in current directory

// Our graph consists of a list of nodes where each node is represented as follows:

double maxChange = 10.00;

class Graph_Node{

private:
    string Node_Name;  // Variable name 
    vector<int> Children; // Children of a particular node - these are index of nodes in graph.
    vector<string> Parents; // Parents of a particular node- note these are names of parents
    int nvalues;  // Number of categories a variable represented by this node can take
    vector<string> values; // Categories of possible values
    vector<double> CPT; // conditional probability table as a 1-d array . Look for BIF format to understand its meaning

public:
    // Constructor- a node is initialised with its name and its categories
    
    vector < int> parentIds;
    vector<int> multiplies;
    vector<int> parentNumValues;
    
    Graph_Node(string name,int n,vector<string> vals){
        Node_Name=name;
        nvalues=n;
        values=vals;
    }


    string get_name(){
        return Node_Name;
    }

    vector<int> get_children(){
        return Children;
    }

    vector<string> get_Parents(){
        return Parents;
    }

    vector<double> get_CPT(){
        return CPT;
    }

    int get_nvalues(){
        return nvalues;
    }

    int get_value_index(string value){
        for(int i=0;i<values.size();i++){
            if(values[i]==value) return i;
        }
    }

    string get_index_value(int n){
        return values[n];
    }

    vector<string> get_values(){
        return values;
    }

    void set_CPT(vector<double> &new_CPT){
        CPT.clear();
        CPT=new_CPT;
    }

    void set_Parents(vector<string> &Parent_Nodes){
        Parents.clear();
        Parents=Parent_Nodes;
    }

    // add another node in a graph as a child of this node
    int add_child(int new_child_index){
        for(int i=0;i<Children.size();i++)
        {
            if(Children[i]==new_child_index)
                return 0;
        }
        Children.push_back(new_child_index);
        return 1;
    }

    double probGivenParent(int k, vector< int > &values){
        // k is the index of value taken by this particular node
        // values is the vector containing indices of values taken by parents and their total number of values possible
        int n = parentIds.size();
        int index = 0;
        for(int i=0;i<n;i++){
            index = index + (values[i])*(multiplies[i+1]);
        }
        index = k*multiplies[0] + index;
        return CPT[index];
    }
};

// The whole network represted as a list of nodes
class network{
    unordered_map<string, int> mapping;

public:
    vector<Graph_Node> Pres_Graph;
    int addNode(Graph_Node node){
        Pres_Graph.push_back(node);
        mapping.insert({node.get_name(), Pres_Graph.size()-1});
        return 0;
    }
    
    void setParentIds(){
        vector<string> parents;
        
        for(int it=0;it<Pres_Graph.size();it++){ // Pres_Graph.size() is the size of allvalues // For this particular model provided
            parents = Pres_Graph[it].get_Parents();
            Pres_Graph[it].parentIds.clear();
            for(int i=0;i<parents.size();i++){
                Pres_Graph[it].parentIds.push_back(mapping[parents[i]]);
            }
        }
    }

    void setParentNumValues(){
        for(int it=0;it<Pres_Graph.size();it++){ // Pres_Graph.size() is the size of allvalues // For this particular model provided
            for(int i=0;i<Pres_Graph[it].parentIds.size();i++){
                int temp = Pres_Graph[Pres_Graph[it].parentIds[i]].get_nvalues();
                Pres_Graph[it].parentNumValues.push_back(temp);
            }
        }
    }

    void setMultiplies(){
        for(int it=0;it<Pres_Graph.size();it++){ // Pres_Graph.size() is the size of allvalues // For this particular model provided
            Pres_Graph[it].multiplies.clear();
            int n = Pres_Graph[it].parentIds.size();
            int temp = 1;
            for(int i=0;i<n;i++){
                temp = temp*(Pres_Graph[it].parentNumValues[i]);
            }
            Pres_Graph[it].multiplies.push_back(temp);
            for(int i=0;i<n;i++){
                temp = temp/(Pres_Graph[it].parentNumValues[i]);
                Pres_Graph[it].multiplies.push_back(temp);
            }
        }
    }

    int netSize(){
        return Pres_Graph.size();
    }

    // get the index of node with a given name
    int get_index(string val_name){
        vector<Graph_Node>::iterator listIt;
        int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return count;
            count++;
        }
        return -1;
    }
    
    // get the node at nth index
    vector<Graph_Node>::iterator get_nth_node(int n){
       vector<Graph_Node>::iterator listIt;
        int count=0;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(count==n)
                return listIt;
            count++;
        }
        return listIt; 
    }
    
    //get the iterator of a node with a given name
    vector<Graph_Node>::iterator search_node(string val_name){
        vector<Graph_Node>::iterator listIt;
        for(listIt=Pres_Graph.begin();listIt!=Pres_Graph.end();listIt++)
        {
            if(listIt->get_name().compare(val_name)==0)
                return listIt;
        }
    
            cout<<"node not found\n";
        return listIt;
    }

    void initializeExpec(vector < vector<double> > &expec, vector< vector<int> > &records, vector<bool> &isMissing){
        int rows = records.size();
        int columns = records[0].size();
        vector<double> v1;
        for(int i=0;i<rows;i++){
            v1.clear();
            if(!isMissing[i]){
                v1.push_back(-1);
                continue;
            }
            for(int j=0;j<columns;j++){
                int a,b;
                a= records[i][j];
                if(a!=-1) continue;
                b = Pres_Graph[j].get_nvalues();
                
                for(int d=0;d<b;d++){
                    v1.push_back(0);
                }
            }
            expec.push_back(v1);
        }
    }

    void findExpectations(int r, int c, vector < vector<double> > &expec, vector< vector<int> > &records){
        // expectation of r,c has to be replaced with new vector
        vector<double> newExpectations;
        int n = Pres_Graph[c].get_nvalues(); // size of the newExpectations
        
        for(int i=0;i<n;i++){
            records[r][c] = i;
            newExpectations.push_back(joint_distribution(records[r]));
            records[r][c] = -1;   
        }
        double sum = 0.0;
        for(int i=0;i<n;i++){
            sum += newExpectations[i];
        }
        for(int i=0;i<n;i++){
            newExpectations[i] /= sum;
        }
        expec[r] = newExpectations;
    }

    void maximize(vector < vector<double> > &expec, vector< vector<int> > &records){

        vector<double> v;
        for(int current=0;current<Pres_Graph.size();current++){
            
            // Update CPT for each node, this loop corresponds to iterating over each node

            v = Pres_Graph[current].get_CPT(); 
            
            int numVars = v.size(); // size of CPT
            
            int numValues = Pres_Graph[current].get_nvalues(); // number of values current value can take

            vector<int> multiplies;
            // int n = parents.size();
            double ttl = Pres_Graph[current].multiplies[0];
            
            vector<double> sumTemp1(numVars, 0);
            vector<double> sumTemp2(numVars, 0);

            int rows = expec.size();
            int columns = expec[0].size();

            int missingId1=-1, missingId2=-1;
            int missingIndex1 = -1, missingIndex2 = -1;
            vector<int> temp1;
            int currentTemp = -1;
            for(int i=0;i<rows;i++){
                temp1.clear();
                missingId1 = -1;
                missingId2 = -1;
                currentTemp = -1;
                missingIndex1 = missingIndex2 = -1;
                int a = records[i][current];
                if(a==-1){
                    missingId1 = 0;
                    missingIndex1 = current;
                }
                currentTemp = a;
                for(int j=0;j<Pres_Graph[current].parentIds.size();j++){
                    a = records[i][Pres_Graph[current].parentIds[j]];
                    if(a==-1){
                        missingIndex1 = missingIndex2 = Pres_Graph[current].parentIds[j];
                        missingId1 =j+1;
                        missingId2 = j;
                    }
                    temp1.push_back(a);
                }
                // now I have got the necessary values from this record to work on

                if(missingId1==-1){
                    int n = Pres_Graph[current].parentIds.size();
                    int index = 0;
                    for(int p=0;p<n;p++){
                        index = index + (temp1[p])*(Pres_Graph[current].multiplies[p+1]);
                    }
                    index = currentTemp*Pres_Graph[current].multiplies[0] + index;
                    sumTemp1[index] +=1;
                }
                else{
                    int n = Pres_Graph[current].parentIds.size();
                    
                    int nval = Pres_Graph[missingIndex1].get_nvalues();
                    double inc = 0.0;
                    for(int d=0;d<nval;d++){
                        inc = expec[i][d];
                        int index = 0;
                        if(missingId1==0){
                            currentTemp = d;
                        }
                        else{
                            temp1[missingId1-1] = d;
                        }

                        for(int p=0;p<n;p++){
                            index = index + (temp1[p])*(Pres_Graph[current].multiplies[p+1]);
                        }
                        index = currentTemp*Pres_Graph[current].multiplies[0] + index;
                        sumTemp1[index] += inc;
                        if(missingId1==0){
                            currentTemp = -1;
                        }
                        else{
                            temp1[missingId1-1] = -1;
                        }
                    } 
                }

                if(missingId2==-1){
                    int n = Pres_Graph[current].parentIds.size();
                    int index = 0;
                    for(int p=0;p<n;p++){
                        index = index + (temp1[p])*(Pres_Graph[current].multiplies[p+1]);
                    }
                    for(int k=0;k<numValues;k++){
                        int index1 = k*Pres_Graph[current].multiplies[0] + index;
                        sumTemp2[index1] +=1;
                    }
                }
                else{
                    int n = Pres_Graph[current].parentIds.size();
                    
                    int nval = Pres_Graph[missingIndex2].get_nvalues();
                    double inc = 0.0;
                    for(int d=0;d<nval;d++){
                        inc = expec[i][d];
                        int index = 0;
                        
                        temp1[missingId2] = d;

                        for(int p=0;p<n;p++){
                            index = index + (temp1[p])*(Pres_Graph[current].multiplies[p+1]);
                        }

                        for(int k=0;k<numValues;k++){
                            int index1 = k*Pres_Graph[current].multiplies[0] + index;
                            sumTemp2[index1] +=inc;
                        }
                        
                        temp1[missingId2] = -1;
                    }
                }
            }

            for(int i=0;i<numVars;i++){
                double db = (sumTemp1[i] + 0.001)/(sumTemp2[i] + 0.001*numValues);
                maxChange = max(maxChange, fabs(v[i]-db));
                v[i] = db;
                if(v[i]<0.0001) v[i] = 0.0001;
                if(v[i]>0.9999) v[i] = 0.9999;
            }

            Pres_Graph[current].set_CPT(v);
        }
    }

    void initializeRandom(){
        srand(time(NULL));
        vector<Graph_Node> :: iterator itr = Pres_Graph.begin();
        while(itr!=Pres_Graph.end()){
            vector<double> cpt = itr->get_CPT();
            int n = cpt.size();
            for(int i=0;i<n;i++){
                if(cpt[i]==-1){
                    cpt[i] = ((rand()%9999)+1)/10000.0; // it won't take 0 and 1
                    
                }
            }
            itr->set_CPT(cpt);
            itr++;
        }
    }

    double joint_distribution(vector<int> &allvalues){
        double result = 1.0;
        vector<int> parentValues;
        vector< pair<int, int> > parentValueIds;
        
        for(int it=0;it<Pres_Graph.size();it++){ // Pres_Graph.size() is the size of allvalues // For this particular model provided
            parentValues.clear();
            
            for(int i=0;i<Pres_Graph[it].parentIds.size();i++){
                parentValues.push_back(allvalues[Pres_Graph[it].parentIds[i]]);
            }
            
            int k = allvalues[it];
            result = result*Pres_Graph[it].probGivenParent(k, parentValues);
        }
        return result;
    }

    void finalPrint(string fileArg){
        ofstream outfile("solved_alarm.bif");
        string temp;
        string line;
        if(outfile.is_open()){
            outfile << fixed;
            outfile.precision(4);
            ifstream myfile(&fileArg[0]);
            if(myfile.is_open()){
                while(!myfile.eof()){
                    stringstream ss;
                    getline(myfile, line);
                    ss.str(line);
                    ss >> temp;
                    
                    if(temp.compare("probability")==0){
                        ss >> temp;
                        ss >> temp; // now temp has the corresponding variable
                        int temp1 = mapping[temp]; // temp1 has the corresponding index
                        outfile << line << endl;
                        getline(myfile, line);
                        outfile << "\ttable";
                        vector<double> v = Pres_Graph[temp1].get_CPT();
                        for(int i=0;i<v.size();i++){
                            if(!myfile.eof()) outfile << " " << v[i];
                        }
                        if(!myfile.eof()) outfile << " ;" << endl;
                    }
                    else{
                        if(!myfile.eof()) outfile << line << endl;
                    }
                }
                myfile.close();
            } 
            else{
                cerr << "alarm.bif couldn't be opened" << endl;
                exit(1);
            }
            outfile.close();
        }
        else{
            cerr<< "output file couldn't be generated" << endl;
            exit(1);
        }
    }
};

network read_network(string fileArg){
    network Alarm;
    string line;
    int find=0;
    ifstream myfile(&fileArg[0]); 
    string temp;
    string name;
    vector<string> values;
    
    if (myfile.is_open())
    {
        while (! myfile.eof() )
        {
            stringstream ss;
            getline (myfile,line);
            
            ss.str(line);
            ss>>temp;
            
            if(temp.compare("variable")==0)
            {
                ss>>name;
                getline (myfile,line);
               
                stringstream ss2;
                ss2.str(line);
                for(int i=0;i<4;i++){
                    ss2>>temp;
                }
                values.clear();
                while(temp.compare("};")!=0)
                {
                    values.push_back(temp);
                    
                    ss2>>temp;
                }
                Graph_Node new_node(name,values.size(),values);
                int pos=Alarm.addNode(new_node);
                    
            }
            else if(temp.compare("probability")==0)
            {    
                ss>>temp;
                ss>>temp;
                
                vector<Graph_Node>::iterator listIt;
                vector<Graph_Node>::iterator listIt1;
                listIt=Alarm.search_node(temp);
                int index=Alarm.get_index(temp);
                ss>>temp;
                values.clear();
                while(temp.compare(")")!=0)
                {
                    listIt1=Alarm.search_node(temp);
                    listIt1->add_child(index);
                    values.push_back(temp);
                    ss>>temp;
                }
                listIt->set_Parents(values);
                getline (myfile,line);
                stringstream ss2;
                
                ss2.str(line);
                ss2>> temp;
                
                ss2>> temp;
                
                vector<double> curr_CPT;
                string::size_type sz;
                while(temp.compare(";")!=0)
                {
                    curr_CPT.push_back(atof(temp.c_str()));
                    ss2>>temp;
                }
                listIt->set_CPT(curr_CPT);
            }   
        }
        if(find==1)
        myfile.close();
    }
    return Alarm;
}


int main(int argc, char *argv[]){
	time_t tinit;
	time(&tinit);
    network Alarm;
    Alarm=read_network(argv[1]);
    ifstream inFile;

    inFile.open(argv[2]);
    if(!inFile){
        cerr<< "File coundn't be opened" << endl;
        exit(1);
    }
    vector< vector<int> > records;
    vector< pair<int, int> > missing;
    vector<bool> isMissing;
    string line;
    int p= 0;
    while(!(inFile.eof())){
        getline(inFile, line);
        stringstream ss;
        ss.str(line);
        string temp;
        vector<int> v;
        int c = 0;
        bool miss = false;
        while(c<Alarm.netSize()){
            ss>>temp;
            
            if(temp=="\"?\""){
                missing.push_back({p, c});
                v.push_back(-1);
                miss = true;
            }
            else v.push_back(Alarm.Pres_Graph[c].get_value_index(temp));
            c++;
        }
        
        records.push_back(v);
        isMissing.push_back(miss);
        p++;
    }
    inFile.close();

    Alarm.setParentIds();
    Alarm.setParentNumValues();
    Alarm.setMultiplies();
    
    // p is the total number of records
    vector< vector<double> > expec;

    // initialize random probabilities 
    Alarm.initializeRandom();

    Alarm.initializeExpec(expec, records, isMissing);

    int count  = 0;
    time_t tstart;
    time(&tstart);

    int restTime = 570 - difftime(tstart, tinit);

    while(maxChange>0.0001){
        // Expectation step
        time_t ttemp;
        time(&ttemp);
        if(difftime(ttemp, tstart)> restTime){
        	break;
        }
        maxChange = 0;
        for(int i=0;i<missing.size();i++){
            // find expectation of all the missing
            int row,column;
            row = missing[i].first; // record number
            column = missing[i].second; // missing attribute
            Alarm.findExpectations(row, column, expec, records);
        }
        
        //Maximization step ( update CPT )
        Alarm.maximize(expec, records);
        
        count++;
    }
    time_t tfinish;
    time(&tfinish);

    cout << "Total running time of the code is " << difftime(tfinish,tinit) << " seconds" << endl;

    Alarm.finalPrint(argv[1]);
    
    cout<<"Perfect! Hurrah! \n";    
}