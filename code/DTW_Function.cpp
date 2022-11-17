#include "DTW_Function.h"

DTW_Function::DTW_Function(int mfcc_num){
    MFCC_P=mfcc_num;
    MFCC_Pf=float(mfcc_num);
}

DTW_Function::~DTW_Function(){
    return;
}

float DTW_Function::ComputeDTW(vector<vector<float> > cep1, vector<vector<float> > cep2)
{
    vector<float> temp;
    for(int i=0;i<cep1.size();i++)
        for(int j=0;j<cep1[i].size();j++)
            temp.push_back(cep1[i][j]);
    int stdlength=temp.size();
    float * stdmfcc = new float[stdlength];
    std::copy(temp.begin(),temp.end(),stdmfcc);

    vector<float> temp1;
    for(int i=0;i<cep2.size();i++)
        for(int j=0;j<cep2[i].size();j++)
            temp1.push_back(cep2[i][j]);
    int testlen=temp1.size();
    float * testmfcc = new float[testlen];
    std::copy(temp1.begin(),temp1.end(),testmfcc);
    return ComputeDTW(stdmfcc,testmfcc,stdlength/MFCC_P,testlen/MFCC_P);
}

float DTW_Function::ComputeDTW(float *cep1, float *cep2, int num1, int num2){
    struct record
    {		int x;
                int y;
    };
    struct point
    {		int x,y;
                float minvalue;
                        int stepnum;
                                bool recheck;
    };
    record * re;
    record * newre;

    newre=new record[num1*num2];
    re=new record[num1*num2];
    int renum;
    int newrenum=0;
    int i,j;
    point * poi;
    poi=new point[num1*num2];

    for(i=0;i<num1*num2;i++)
    {
        poi[i].recheck=0;
        poi[i].minvalue=INF;
        poi[i].stepnum=0;
    }

    for(i=0;i<5;i++)
    {
        if(i==0)  {	re[i].x=1; re[i].y=1; }
        if(i==1  && (num1 * num2 > 1) )  {	re[i].x=1; re[i].y=2; }
        if(i==2  && (num1 * num2 > 2) )  {	re[i].x=1; re[i].y=3; }
        if(i==3  && (num1 * num2 > 3) )  {	re[i].x=2; re[i].y=1; }
        if(i==4  && (num1 * num2 > 4) )  {	re[i].x=3; re[i].y=1; }
        if(num1 * num2 > 5){
        poi[(re[i].y-1)*num1+re[i].x-1].minvalue=Distance(cep1,cep2,re[i].x,re[i].y);
        poi[(re[i].y-1)*num1+re[i].x-1].stepnum=1;
        }
    }
    renum=5;
    int newx,newy;                   //newvalue;
    for(i=0;i<renum;i++)
    {
        for(j=0;j<3;j++)
        {
            if(j==0){ newx=re[i].x+1; newy=re[i].y+2; }
            if(j==1){ newx=re[i].x+1; newy=re[i].y+1; }
            if(j==2){ newx=re[i].x+2; newy=re[i].y+1; }

            if(newx>=num1||newy>=num2)
                continue;
            if(fabs(newx-newy)<=fabs(num1-num2)+3)
            {
                if(poi[(newy-1)*num1+newx-1].recheck==0)
                {
                    newre[newrenum].x=newx;
                    newre[newrenum].y=newy;
                    newrenum++;
                }
                float tmpdis;
                int addstepnum;
                if(j==0){ tmpdis=Distance(cep1,cep2,newx-1,newy-1)*2+Distance(cep1,cep2,newx,newy); addstepnum=2;}
                if(j==1){ tmpdis=Distance(cep1,cep2,newx,newy)*2; addstepnum=1;}
                if(j==2){ tmpdis=Distance(cep1,cep2,newx-1,newy-1)*2+Distance(cep1,cep2,newx,newy); addstepnum=2;}
                if(poi[(newy-1)*num1+newx-1].minvalue>(poi[(re[i].y-1)*num1+re[i].x-1].minvalue+tmpdis))
                {
                    poi[(newy-1)*num1+newx-1].minvalue=(poi[(re[i].y-1)*num1+re[i].x-1].minvalue+tmpdis);
                    poi[(newy-1)*num1+newx-1].stepnum=poi[(re[i].y-1)*num1+re[i].x-1].stepnum+addstepnum;
                }
                if(poi[(newy-1)*num1+newx-1].recheck==0)
                    poi[(newy-1)*num1+newx-1].recheck=1;
            }
        }
        if(newrenum!=0 && i>=(renum-1))
        {
            renum=newrenum;
            newrenum=0;
            struct	record * tt;
            tt=re;
            re=newre;
            newre=tt;
            i=-1;
        }
    }
    float min=INF;
    for(j=0;j<renum;j++)
    {
        if((poi[(re[j].y-1)*num1+re[j].x-1].minvalue)/poi[(re[j].y-1)*num1+re[j].x-1].stepnum<min)
            min=(poi[(re[j].y-1)*num1+re[j].x-1].minvalue)/poi[(re[j].y-1)*num1+re[j].x-1].stepnum;
    }

    //	min;
    delete []poi;
    delete []newre;
    delete []re;
    delete []cep1;
    delete []cep2;

    return min;
}

float DTW_Function::Distance(float *ps1, float *ps2, int k1, int k2){
    int i=0;
    float sum=0;
    for(i=0;i<MFCC_P;i++)
        sum+=(1+MFCC_Pf/2*(float)sin(PI*i/MFCC_Pf))*(ps1[k1+i]-ps2[k2+i])*(ps1[k1+i]-ps2[k2+i]);

    return sum;
}

//MY FUNCTION

vector<string> DTW_Function::split_str(string str)
{
    string word = "";
    vector<string> temp;
    for (auto x : str)
    {
        if (x == ' ')
        {
            temp.push_back(word);
            word = "";
        }
        else {
            word = word + x;
        }
    }
    return temp;
}

vector<vector<string> >   DTW_Function::read_mfcc_from_csv(string filename){
    ifstream  data(filename);
    string line;
    vector< string > readcsv;
    vector< string > line_csv_temp;
    vector<vector<string> > temp;
    temp.clear();
    //vector<vector<vector<float> > > mfcc_feat;
    while(getline(data,line))
    {
        stringstream lineStream(line);
        string cell;
        vector<string> parsedRow;
        while(getline(lineStream,cell,','))
        {
            parsedRow.push_back(cell);
        }
        parsedRow[3].push_back(' ');
        readcsv.push_back(parsedRow[3]);
    }
    readcsv.erase(readcsv.begin());


    for (int i = 0; i < readcsv.size(); i++)    //total samples number
    {
        line_csv_temp = split_str(readcsv[i]);      //MFCC feature flatten
        temp.push_back(line_csv_temp);
    }
    return temp;
}

vector<vector<float> > DTW_Function::rebuilt_mfcc_feat(vector<string> flatten_mfcc){  //turn a flatten string to a 13*frame_num vector
     vector<vector<float> > temp;
     vector<float> col_temp;
     int row_tag = 0;
     int col_tag = 0;
     for (int j = 0; j < (flatten_mfcc.size() / 13); j++)
     {
        for (int i = 0; i < 13; i++)
        {
            col_temp.push_back(stof(flatten_mfcc[j*13 + i]));
           //row_tag = i % 13;
           //col_tag = int(i / 13);
           //temp[row_tag][col_tag] = /* stof(flatten_mfcc[i]) */0;
        }
        temp.push_back(col_temp);
        col_temp.clear();
     }
     return temp;
}

int DTW_Function::iteration_dtw(vector<vector<string> > match_csv_content, vector<vector<float> >  query_mfcc_feat, float* min_dist, int* min_idx){
    float normal_dist;
    vector<float> eud_distance;
    vector<vector<float> >  match_mfcc_feat;
    for (int i = 0; i < match_csv_content.size(); i++)
    {
        if((match_csv_content[i].size() % 13) != 0){
            cout << "WRONG READ!" << "MFCC Length is: "<< match_csv_content[i].size() << endl;
            /* for (int j = 0; j < match_csv_content[i].size(); j++)
            {
                cout << match_csv_content[i][j] << endl;
            } */
            break;
        }
        else{
            float distance = 0;
            float normal_distance = 0;
            match_mfcc_feat = rebuilt_mfcc_feat(match_csv_content[i]);
            distance = ComputeDTW(query_mfcc_feat,match_mfcc_feat);
            normal_dist = distance / (query_mfcc_feat.size() + match_mfcc_feat.size());
            eud_distance.push_back(normal_dist);
        }
    }
    //eud_distance = standard_normal(eud_distance);   //no standard
    auto    min_pos   = std::min_element(eud_distance.begin(), eud_distance.end());
    *min_idx = (min_pos - eud_distance.begin());
    *min_dist = *min_pos;

    return 0;
}

vector<float> DTW_Function::standard_normal(vector<float> ini_vector){
    vector<float> temp;
    float mean;
    float std_var;
    mean = vector_mean(ini_vector);
    std_var = vector_std_variance(ini_vector, mean);
    for (int i = 0; i < ini_vector.size(); i++)
    {
        temp.push_back((ini_vector[i] - mean) / std_var);
    }
    return temp;
}

float DTW_Function::vector_mean(vector<float> distance){
    return std::accumulate(distance.begin(), distance.end(), 0) / distance.size();
}

float DTW_Function::vector_std_variance(vector<float> distance, float mean){
    float temp = 0;
    for (int i = 0; i < distance.size(); i++)
    {
        temp += pow((distance[i] - mean), 2);
    }
    temp = sqrt(temp / distance.size());
    return temp;
}