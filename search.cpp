// Do NOT add any other includes
#include "search.h"

SearchEngine::SearchEngine(){
    // Implement your function here 

}

SearchEngine::~SearchEngine(){
    // Implement your function here  
    for(int i=0;i<sentence_info.size();i++){
        delete sentence_info[i];
    }
    
}


void SearchEngine::insert_sentence(int book_code, int page, int paragraph, int sentence_no, string sentence){
    // Implement your function here 
    string lower_cased_sentence;
    for(int i=0;i<sentence.length();i++){
        if(int(sentence[i]>64 && int(sentence[i])<91)){
            lower_cased_sentence+=char(int(sentence[i]+32));
        }
        else{
            lower_cased_sentence+=sentence[i];
        }
    }
    text_data.push_back(lower_cased_sentence);
    sentence_info.push_back(new Node(book_code,page,paragraph,sentence_no,0));
    return;
}
vector<int> KMP(string input,string data,vector<int> &v){
    vector<int> result;
    if(input.length()>data.length()){
        return result;
    }
    int i = 0; // index for txt[]
    int j = 0; // index for pat[]
    int n=data.length();
    int m=input.length();
    while ((n-i) >= (m- j)) {
        if (input[j] == data[i]) {
            j++;
            i++;
        }
        if (j == m) {
            result.push_back(i-j);
            j = v[j - 1];
        }
        else if (i < n && input[j] != data[i]) {
            if (j != 0)
                j = v[j - 1];
            else
                i = i + 1;
        }
    }

    return result;
}
Node* SearchEngine::search(string pattern, int& n_matches){
   
    if(pattern.length()==0){
        return nullptr;
    }
    //cout<<"here"<<endl;
    int start=text_data.size();
    n_matches=0;
    vector<int> v(pattern.length());
    int len = 0;
    v[0] = 0;
    int i = 1;
    
    string lower_pattern="";
    for(int i=0;i<pattern.length();i++){
        if(int(pattern[i])>64 && int(pattern[i])<91){
            lower_pattern+=char(int(pattern[i]+32));
        }
        else{
            lower_pattern+=pattern[i];
        }
    }
    pattern=lower_pattern;
    while (i < pattern.length()) {
        if (pattern[i] == v[len]) {
            len++;
            v[i] = len;
            i++;
        }
        else{
            if (len != 0){
                len = v[len - 1];
            }
            else{
                v[i] = 0;
                i++;
            }
        }
    }
    
    vector<Node *> temp;
    
    for(int i=0;i<text_data.size();i++){
        vector<int> res=KMP(pattern,text_data[i],v);
        for(int j=0;j<res.size();j++){
            if(temp.size()==0){
                //cout<<res[j]<<endl;
                
                temp.push_back(new Node(sentence_info[i]->book_code,sentence_info[i]->page,sentence_info[i]->paragraph,sentence_info[i]->sentence_no,res[j]));
                head_of_search.push_back(temp[0]);
                temp[n_matches]->left=nullptr;
                temp[n_matches]->right=nullptr;
                n_matches++;
            }
            else{
                temp.push_back(new Node(sentence_info[i]->book_code,sentence_info[i]->page,sentence_info[i]->paragraph,sentence_info[i]->sentence_no,res[j]));
                temp[n_matches]->left=nullptr;
                temp[n_matches]->left=temp[n_matches-1];
                temp[n_matches-1]->right=temp[n_matches];
                n_matches++;

            }
            

            
        }
        
    }
    
    //cout<<n_matches<<endl;
    if(n_matches==0){
        return nullptr;
    }
    
    //cout<<n_matches<<endl;
    temp[0]->left=nullptr;
    temp[n_matches-1]->right=nullptr;
    return temp[0];
}
// #include <bits/stdc++.h>
// #include <chrono>
// using namespace std;

// int main(){

//     srand(time(0));
   
//     SearchEngine* E=new SearchEngine();
//     vector<char> characters;
//     for (char c = 'a'; c <= 'z'; ++c) {
//         characters.push_back(c);
//     }

//     for (char c = 'A'; c <= 'Z'; ++c) {
//         characters.push_back(c);
//     }
//     characters.push_back(';');
//     characters.push_back(',');
//     vector<string> m;
    
//     for(int i=0;i<1000;i++){
//         string sentence="";
        
//         for(int i=0;i<200;i++){
//             sentence+=char(characters[rand()%characters.size()]);
//             if(rand()%4==0){
//                 sentence+=" ";
//             }
//         }
//         E->insert_sentence(rand()%100,rand()%100,rand()%100,rand(),sentence);
//         for(int i=0;i<3;i++){
//             int p=rand()%sentence.length();
//             int q=rand()%sentence.length();
//             m.push_back(sentence.substr(min(p,q),max(p,q)));
//         }
//         for(int i=0;i<3;i++){
//             string serachable_string="";
//             for(int j=0;j<10;j++){
//                 serachable_string+=char(characters[rand()%characters.size()]);
//                 if(rand()%4==0){
//                     serachable_string+="";
//                 }
//             }
//             m.push_back(serachable_string);
//         }
//         //cout<<sentence <<endl;
//     }
//     cout<<"inserted all"<<endl;

//     cout<<m.size()<<endl;
    
//     for(int i=0;i<m.size();i++){
//         //cout<<"|"<<flush;
//         int r=0;
//         int g=0;
//         Node* p=E->search(m[i],r);
        
        
//         while(p!=nullptr){
//             g++;
//             //cout<<p->offset<<endl;
//             p=p->right;
            
//         }
        
        
//         if(r!=g){
            
//             std::cout<<"Fishy fish is not good fish "<<r<<" "<<g<<endl;
//             cout<<m[i]<<" "<<m[i].length()<<endl;
            
//         }
        
        
//     }
//     delete E;

    
//     return 0;

// }
//*/
/*
int main(){
    SearchEngine* E=new SearchEngine();
    E->insert_sentence(0,0,0,0,"abcdefghijkabcdefghijk");
    E->insert_sentence(0,0,0,0,"abcdef");
    E->insert_sentence(0,0,0,0,"a");
    E->insert_sentence(0,0,0,0,"aaaaaaaaaa");
    for(int i=0;i<100;i++){
        string f="";
        f+="xe";
        E->insert_sentence(0,0,0,0,f);
    }

    int r=0;
    int g=0;
    E->search("abcde",r);
    E->search("x",g);
    cout<<r<<endl;
    cout<<g<<endl;
    return 0;
}
 //*/
