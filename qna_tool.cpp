#include "qna_tool.h"

#include <assert.h>
#include <algorithm>
#include <sstream>

using namespace std;

string stemmedword(string word) {
    return word;
}

#define FILENAME "mahatma-gandhi-collected-works-volume-1.txt"

pair<Trie *, Trie *> parse_csv(string csv_file) {
    Trie *t = new Trie();

    Trie *stemt = new Trie();  // stemmed trie

    ifstream csv;
    csv.open(csv_file);
    // ignore first line
    string junk;
    csv >> junk;
    while (!csv.eof()) {
        string word;
        long long freq;
        string separator;
        string line;
        csv >> line;
        if (line.empty()) {
            continue;
        }
        int i = 0;
        while (line[i] != ',') {
            word += line[i];
            i++;
        }
        i++;
        freq = stoll(line.substr(i, line.size() - i));

        t->set_freq(word, freq);

        stemt->set_freq(stemmedword(word), freq);  // stemmed trie
    }
    return {t, stemt};
}

QNA_tool::QNA_tool() {
    pair<Trie *, Trie *> p = parse_csv("unigram_freq.csv");
    StandardText = p.first;

    StemmedStandardText = p.second;  // stemmed trie
}

QNA_tool::~QNA_tool() {
    delete StandardText;

    delete StemmedStandardText;  // stemmed trie
}

void heapify(vector<pair<double, Node *>> &t, int n, int i) {
    int k = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    // cout<<"left "<<left<<endl;
    // cout<<"right "<<right<<endl;
    if (left < n && t[left].first > t[k].first) {
        k = left;
    }
    if (right < n && t[right].first > t[k].first) {
        k = right;
    }
    if (k != i) {
        pair<double, Node *> temp = t[i];
        t[i] = t[k];
        t[k] = temp;
        heapify(t, n, k);
    }
}
void buildHeap(vector<pair<double, Node *>> &t, int n) {
    int start = (n / 2) - 1;
    for (int i = start; i > -1; i--) {
        heapify(t, n, i);
    }
}
void heapify_down(vector<pair<double, Node *>> &t) {
    int now = 0;
    int n = t.size();

    while (now < n) {
        // cout<<now<<" "<<flush;
        int left = 2 * now + 1;
        int right = 2 * now + 2;
        if (left < n && right < n) {
            if (t[left].first < t[right].first) {
                if (t[now].first < t[right].first) {
                    pair<double, Node *> temp = t[now];
                    t[now] = t[right];
                    t[right] = temp;
                }
                now = right;

            } else {
                if (t[now].first < t[left].first) {
                    pair<double, Node *> temp = t[now];
                    t[now] = t[left];
                    t[left] = temp;
                }
                now = left;
            }

        }

        else if (left < n) {
            if (t[now].first < t[left].first) {
                pair<double, Node *> temp = t[now];
                t[now] = t[left];
                t[left] = temp;
            }
            now = left;

        }

        else {
            return;
        }
    }
}
Node *get_top(vector<pair<double, Node *>> &t) {
    Node *top = t[0].second;

    t[0] = t[t.size() - 1];
    t.pop_back();

    heapify_down(t);

    return top;
}
Node *k_top(vector<pair<double, Node *>> &t, int k) {
    buildHeap(t, t.size());
    if (k > t.size()) {
        k = t.size();
    }

    // cout<<0<<endl;

    Node *root = get_top(t);
    root->left = nullptr;
    Node *now = root;
    for (int i = 1; i < k; i++) {
        // cout<<i<<endl;

        now->right = get_top(t);
        now->right->left = now;
        now = now->right;
    }

    now->right = nullptr;

    return root;
}

Node *QNA_tool::get_top_k_para(string question, int k) {
    // separating words in question
    vector<string> question_words;
    string word = "";
    for (int i = 0; i < question.size(); i++) {
        if (is_separator(question[i])) {
            if (word != "") {
                word = to_lower(word);  // converting to lower case
                question_words.push_back(word);
            }
            word = "";
        } else {
            word += question[i];
        }
    }
    if (word != "") {
        word = to_lower(word);
        question_words.push_back(word);
    }

    // making a trie of question words with frequency equal to thier frequency
    // in csv and corpus for faster score computation

    Trie *question_scores_csv = new Trie();
    Trie *question_scores_corpus = new Trie();

    for (string s : question_words) {
        question_scores_csv->set_freq(s, StandardText->get_freq(s));
        question_scores_corpus->set_freq(s, corpus->get_freq(s));
    }

    vector<pair<double, Node *>> ParaScores;

    for (long long book_num = 0; book_num < ParaTrie.size(); book_num++) {
        for (long long page = 0; page < ParaTrie[book_num].size(); page++) {
            for (long long para = 0; para < ParaTrie[book_num][page].size();
                 para++) {
                double score = 0;
                Trie *t = ParaTrie[book_num][page][para];
                for (string s : question_words) {
                    double numerator = question_scores_corpus->get_freq(s) + 1;
                    double denominator = question_scores_csv->get_freq(s) + 1;
                    double f = t->get_freq(s);
                    score += f * (numerator / denominator);
                }
                Node *n = new Node();
                n->book_code = book_num + 1;
                n->page = page;
                n->paragraph = para;
                pair<double, Node *> p = {score, n};
                ParaScores.push_back(p);
            }
        }
    }
    // cout << "PRINTING VECTOR:" << endl;
    // for (auto i : ParaScores) {
    //     cout << i.first << " "
    //          << " " << i.second->book_code << " " << i.second->page << " "
    //          << i.second->paragraph << endl;
    // }
    // cout << "------------" << endl;
    Node *head_of_k_maximum =
        new Node();  // a function that will give head of k maximum Nodes
    head_of_k_maximum = k_top(ParaScores, k);
    return head_of_k_maximum;
}
string lower_case(string s){
    string result="";
    for(int i=0;i<s.length();i++){
        if(int(s[i])>64 && int(s[i])<91){
            result+=char(int(s[i])+32);

        }
        else{
            result+=s[i];
        }
    }
    return result;
}
string reduced_querry_string(string s){
    //vector<string> not_used={};
    vector<string> not_used={"?","view","s","s","views", "gandhi", "mahatma", "a", "about", "above", "after", "although", "am", "among", "an", "and", "are", "around", "as", "at", "be", "been", "before", "being", "below", "beside", "between", "beyond", "but", "by", "can", "could", "did", "do", "does", "due", "even", "for", "from", "had", "has", "have", "he", "hence", "how", "however", "i", "if", "in", "instead", "into", "is", "nor", "not", "of", "off", "on", "onto", "or", "other", "over", "rather", "reason", "shall", "she", "should", "so", "such", "than", "that", "the", "their", "them", "then", "there", "therefore", "they", "this", "though", "through", "thus", "till", "to", "toward", "under", "unlike", "until", "was", "well", "were", "what", "when", "where", "whereas", "whether", "which", "while", "who", "whom", "whome", "whose", "why", "will", "with", "without", "would", "yet", "you"};

    vector<string> words;
    string word="";
    for(int i=0;i<s.length();i++){
        //!cout<<word<<endl;
        if(s[i]=='.' || s[i]==','  || s[i]=='-'  || s[i]==':'  || s[i]=='\"'  || s[i]=='\''  || s[i]=='!' || s[i]=='(' || s[i]==')'  || s[i]=='?' || s[i]=='[' || s[i]==']' || s[i]==';' || s[i]=='@' || s[i]== ' '){
            //cout<<word<<endl;

            if(word.size()>0){
                
                int find=true;
                word=lower_case(word);
                //cout<<word<<endl;
                
                
                for(int j=0;j<not_used.size();j++){
                    if(not_used[j]==word){
                        find=false;
                    }
                }
                if(find){
                    words.push_back(word);
                }
                
            }
            word="";
        }
        else{

        
            word+=s[i];
        }
    }
    if(word.size()>0){
                
        int find=true;
        
        for(int j=0;j<not_used.size();j++){
            if(not_used[j]==word){
                find=false;
            }
        }
        if(find){
            words.push_back(word);
        }
        
    }
    
    string result="";
    for(auto i:words){
        result+=i+" ";
    }
    return result;
}

void set_same_values(Node *a,Node *b){
    a->book_code=b->book_code;
    a->offset=b->offset;
    a->sentence_no=b->sentence_no;
    a->page=b->page;
    a->paragraph=b->paragraph;
    return;
}



void QNA_tool::query(string question, string filename) {
    // Implement your function here
    string gptquerry = question;
    question=reduced_querry_string(question);
    // separating words in question
    vector<string> Stemmedquestion_words;
    string word = "";
    for (int i = 0; i < question.size(); i++) {
        if (is_separator(question[i])) {
            if (word != "") {
                word = to_lower(word);
                  // converting to lower case
                if(StemmedCorpus->get_freq(word)>0){
                    Stemmedquestion_words.push_back(word);
                }
                
            }
            word = "";
        } else {
            word += question[i];
        }
    }
    if (word != "") {
        word = to_lower(word);
        Stemmedquestion_words.push_back(word);
    }

    for (int i = 0; i < Stemmedquestion_words.size();
         i++) {  // Stemmedquestion_words contains words of question in stemmed
                 // form
        string word = Stemmedquestion_words[i];
        Stemmedquestion_words[i] = stemmedword(word);
    }

    // making a trie of question words with frequency equal to thier frequency
    // in csv and corpus for faster score computation
    

    if(Stemmedquestion_words.size()==0){
        Stemmedquestion_words.push_back("gandhi");
        

        int k = 7;

        Trie *StemmedQuestion_scores_csv = new Trie();
        Trie *StemmedQuestion_scores_corpus = new Trie();

        for (string s : Stemmedquestion_words) {
            StemmedQuestion_scores_csv->set_freq(s,
                                                StemmedStandardText->get_freq(s));
            StemmedQuestion_scores_corpus->set_freq(s, StemmedCorpus->get_freq(s));
        }

        vector<pair<double, Node *> > ParaScores0;
        for (long long book_num = 0; book_num < ParaTrie.size(); book_num++) {
            for (long long page = 0; page < ParaTrie[book_num].size(); page++) {
                for (long long para = 0; para < ParaTrie[book_num][page].size();
                    para++) {
                    double score = 0;
                    Trie *t = ParaTrie[book_num][page][para];
                    string s="gandhi";
                    double numerator =
                        StemmedQuestion_scores_corpus->get_freq(s) + 1;
                    double denominator =
                        StemmedQuestion_scores_csv->get_freq(s) + 1;
                    double f = t->get_freq(s);
                    score += f * (numerator / denominator);
                    
                    Node *n = new Node();
                    n->book_code = book_num + 1;
                    n->page = page;
                    n->paragraph = para;
                    pair<double, Node *> p = {score, n};
                    ParaScores0.push_back(p);
                }
            }
        }
        Node *head_of_k_maximum =new Node();
        head_of_k_maximum=k_top(ParaScores0,7);

        query_llm(filename, head_of_k_maximum,7, "hello", gptquerry);
    }
    else if(Stemmedquestion_words.size()==1){
        vector<pair<int,string> > word_freq;
        
        for (string s : Stemmedquestion_words) {
            word_freq.push_back(make_pair(StemmedCorpus->get_freq(s),s));
            
        }
        sort(word_freq.begin(),word_freq.end());
        for(auto i:word_freq){
            cout<<i.second<<" "<<i.first<<endl;
        }


        



        int k = 10;

        Trie *StemmedQuestion_scores_csv = new Trie();
        Trie *StemmedQuestion_scores_corpus = new Trie();

        for (string s : Stemmedquestion_words) {
            StemmedQuestion_scores_csv->set_freq(s,
                                                StemmedStandardText->get_freq(s));
            StemmedQuestion_scores_corpus->set_freq(s, StemmedCorpus->get_freq(s));
        }

        vector<pair<double, Node *> > ParaScores1;

        for (long long book_num = 0; book_num < ParaTrie.size(); book_num++) {
            for (long long page = 0; page < ParaTrie[book_num].size(); page++) {
                for (long long para = 0; para < ParaTrie[book_num][page].size();
                    para++) {
                    double score = 0;
                    Trie *t = ParaTrie[book_num][page][para];
                    string s=word_freq[0].second;
                    double numerator =
                        StemmedQuestion_scores_corpus->get_freq(s) + 1;
                    double denominator =
                        StemmedQuestion_scores_csv->get_freq(s) + 1;
                    double f = t->get_freq(s);
                    score += f * (numerator / denominator);
                    
                    Node *n = new Node();
                    n->book_code = book_num + 1;
                    n->page = page;
                    n->paragraph = para;
                    pair<double, Node *> p = {score, n};
                    ParaScores1.push_back(p);
                }
            }
        }
        
        Node *head_of_k_maximum1 =
            new Node();  // a function that will give head of k maximum Nodes
        head_of_k_maximum1 = k_top(ParaScores1, k);
        

        Node *head_of_k_maximum =
            new Node();

        int i=0;
        while(i<7 && head_of_k_maximum1!=nullptr){
            set_same_values(head_of_k_maximum,head_of_k_maximum1);
            head_of_k_maximum->right=new Node();
            head_of_k_maximum->right->left=head_of_k_maximum;
            head_of_k_maximum=head_of_k_maximum->right;
            head_of_k_maximum1=head_of_k_maximum1->right;
            i++;
        }
        
        while(head_of_k_maximum->left!=nullptr){
            head_of_k_maximum=head_of_k_maximum->left;
        }
        /* while(head_of_k_maximum->right!=nullptr){
            cout<<head_of_k_maximum->book_code<<endl;
            head_of_k_maximum=head_of_k_maximum->right;
        } */
        std::cout<<".............."<<question<<endl;
        
        query_llm(filename, head_of_k_maximum, 7, "hello", gptquerry);
    }
    else if(Stemmedquestion_words.size()==2){
            vector<pair<int,string> > word_freq;
        
        for (string s : Stemmedquestion_words) {
            word_freq.push_back(make_pair(StemmedCorpus->get_freq(s),s));
            
        }
        sort(word_freq.begin(),word_freq.end());
        for(auto i:word_freq){
            cout<<i.second<<" "<<i.first<<endl;
        }

        



        int k = 10;

        Trie *StemmedQuestion_scores_csv = new Trie();
        Trie *StemmedQuestion_scores_corpus = new Trie();

        for (string s : Stemmedquestion_words) {
            StemmedQuestion_scores_csv->set_freq(s,
                                                StemmedStandardText->get_freq(s));
            StemmedQuestion_scores_corpus->set_freq(s, StemmedCorpus->get_freq(s));
        }

        vector<pair<double, Node *> > ParaScores1;

        for (long long book_num = 0; book_num < ParaTrie.size(); book_num++) {
            for (long long page = 0; page < ParaTrie[book_num].size(); page++) {
                for (long long para = 0; para < ParaTrie[book_num][page].size();
                    para++) {
                    double score = 0;
                    Trie *t = ParaTrie[book_num][page][para];
                    string s=word_freq[0].second;
                    double numerator =
                        StemmedQuestion_scores_corpus->get_freq(s) + 1;
                    double denominator =
                        StemmedQuestion_scores_csv->get_freq(s) + 1;
                    double f = t->get_freq(s);
                    score += f * (numerator / denominator);
                    
                    Node *n = new Node();
                    n->book_code = book_num + 1;
                    n->page = page;
                    n->paragraph = para;
                    pair<double, Node *> p = {score, n};
                    ParaScores1.push_back(p);
                }
            }
        }
        
        Node *head_of_k_maximum1 =
            new Node();  // a function that will give head of k maximum Nodes
        head_of_k_maximum1 = k_top(ParaScores1, k);
        


        vector<pair<double, Node *> > ParaScores2;

        for (long long book_num = 0; book_num < ParaTrie.size(); book_num++) {
            for (long long page = 0; page < ParaTrie[book_num].size(); page++) {
                for (long long para = 0; para < ParaTrie[book_num][page].size();
                    para++) {
                    double score = 0;
                    Trie *t = ParaTrie[book_num][page][para];
                    string s=word_freq[1].second;
                    double numerator =
                        StemmedQuestion_scores_corpus->get_freq(s) + 1;
                    double denominator =
                        StemmedQuestion_scores_csv->get_freq(s) + 1;
                    double f = t->get_freq(s);
                    score += f * (numerator / denominator);
                    
                    Node *n = new Node();
                    n->book_code = book_num + 1;
                    n->page = page;
                    n->paragraph = para;
                    pair<double, Node *> p = {score, n};
                    ParaScores2.push_back(p);
                }
            }
        }
        
        Node *head_of_k_maximum2 =
            new Node();  // a function that will give head of k maximum Nodes
        head_of_k_maximum2 = k_top(ParaScores2, k);

        

        


        Node *head_of_k_maximum =
            new Node();

        int i=0;
        while(i<5 && head_of_k_maximum1!=nullptr){
            set_same_values(head_of_k_maximum,head_of_k_maximum1);
            head_of_k_maximum->right=new Node();
            head_of_k_maximum->right->left=head_of_k_maximum;
            head_of_k_maximum=head_of_k_maximum->right;
            head_of_k_maximum1=head_of_k_maximum1->right;
            i++;
        }
        while(i<7 && head_of_k_maximum2!=nullptr ){
            set_same_values(head_of_k_maximum,head_of_k_maximum2);
            head_of_k_maximum->right=new Node();
            head_of_k_maximum->right->left=head_of_k_maximum;
            head_of_k_maximum=head_of_k_maximum->right;
            head_of_k_maximum2=head_of_k_maximum2->right;
            i++;
        }
        
        while(head_of_k_maximum->left!=nullptr){
            head_of_k_maximum=head_of_k_maximum->left;
        }
        /* while(head_of_k_maximum->right!=nullptr){
            cout<<head_of_k_maximum->book_code<<endl;
            head_of_k_maximum=head_of_k_maximum->right;
        } */
        std::cout<<".............."<<question<<endl;
        
        query_llm(filename, head_of_k_maximum, 7, "hello", gptquerry);
    }
    else{

        vector<pair<int,string> > word_freq;
        
        for (string s : Stemmedquestion_words) {
            word_freq.push_back(make_pair(StemmedCorpus->get_freq(s),s));
            
        }
        sort(word_freq.begin(),word_freq.end());
        for(auto i:word_freq){
            cout<<i.second<<" "<<i.first<<endl;
        }

        



        int k = 10;

        Trie *StemmedQuestion_scores_csv = new Trie();
        Trie *StemmedQuestion_scores_corpus = new Trie();

        for (string s : Stemmedquestion_words) {
            StemmedQuestion_scores_csv->set_freq(s,
                                                StemmedStandardText->get_freq(s));
            StemmedQuestion_scores_corpus->set_freq(s, StemmedCorpus->get_freq(s));
        }

        vector<pair<double, Node *> > ParaScores1;

        for (long long book_num = 0; book_num < ParaTrie.size(); book_num++) {
            for (long long page = 0; page < ParaTrie[book_num].size(); page++) {
                for (long long para = 0; para < ParaTrie[book_num][page].size();
                    para++) {
                    double score = 0;
                    Trie *t = ParaTrie[book_num][page][para];
                    string s=word_freq[0].second;
                    double numerator =
                        StemmedQuestion_scores_corpus->get_freq(s) + 1;
                    double denominator =
                        StemmedQuestion_scores_csv->get_freq(s) + 1;
                    double f = t->get_freq(s);
                    score += f * (numerator / denominator);
                    
                    Node *n = new Node();
                    n->book_code = book_num + 1;
                    n->page = page;
                    n->paragraph = para;
                    pair<double, Node *> p = {score, n};
                    ParaScores1.push_back(p);
                }
            }
        }
        
        Node *head_of_k_maximum1 =
            new Node();  // a function that will give head of k maximum Nodes
        head_of_k_maximum1 = k_top(ParaScores1, k);
        


        vector<pair<double, Node *> > ParaScores2;

        for (long long book_num = 0; book_num < ParaTrie.size(); book_num++) {
            for (long long page = 0; page < ParaTrie[book_num].size(); page++) {
                for (long long para = 0; para < ParaTrie[book_num][page].size();
                    para++) {
                    double score = 0;
                    Trie *t = ParaTrie[book_num][page][para];
                    string s=word_freq[1].second;
                    double numerator =
                        StemmedQuestion_scores_corpus->get_freq(s) + 1;
                    double denominator =
                        StemmedQuestion_scores_csv->get_freq(s) + 1;
                    double f = t->get_freq(s);
                    score += f * (numerator / denominator);
                    
                    Node *n = new Node();
                    n->book_code = book_num + 1;
                    n->page = page;
                    n->paragraph = para;
                    pair<double, Node *> p = {score, n};
                    ParaScores2.push_back(p);
                }
            }
        }
        
        Node *head_of_k_maximum2 =
            new Node();  // a function that will give head of k maximum Nodes
        head_of_k_maximum2 = k_top(ParaScores2, k);

        vector<pair<double, Node *> > ParaScores3;

        for (long long book_num = 0; book_num < ParaTrie.size(); book_num++) {
            for (long long page = 0; page < ParaTrie[book_num].size(); page++) {
                for (long long para = 0; para < ParaTrie[book_num][page].size();
                    para++) {
                    double score = 0;
                    Trie *t = ParaTrie[book_num][page][para];
                    string s=word_freq[2].second;
                    double numerator =
                        StemmedQuestion_scores_corpus->get_freq(s) + 1;
                    double denominator =
                        StemmedQuestion_scores_csv->get_freq(s) + 1;
                    double f = t->get_freq(s);
                    score += f * (numerator / denominator);
                    
                    Node *n = new Node();
                    n->book_code = book_num + 1;
                    n->page = page;
                    n->paragraph = para;
                    pair<double, Node *> p = {score, n};
                    ParaScores3.push_back(p);
                }
            }
        }
        
        Node *head_of_k_maximum3 =
            new Node();  // a function that will give head of k maximum Nodes
        head_of_k_maximum3 = k_top(ParaScores3, k);


        Node *head_of_k_maximum =
            new Node();

        int i=0;
        while(i<4 && head_of_k_maximum1!=nullptr){
            set_same_values(head_of_k_maximum,head_of_k_maximum1);
            head_of_k_maximum->right=new Node();
            head_of_k_maximum->right->left=head_of_k_maximum;
            head_of_k_maximum=head_of_k_maximum->right;
            head_of_k_maximum1=head_of_k_maximum1->right;
            i++;
        }
        while(i<6 && head_of_k_maximum2!=nullptr ){
            set_same_values(head_of_k_maximum,head_of_k_maximum2);
            head_of_k_maximum->right=new Node();
            head_of_k_maximum->right->left=head_of_k_maximum;
            head_of_k_maximum=head_of_k_maximum->right;
            head_of_k_maximum2=head_of_k_maximum2->right;
            i++;
        }
        while(i<7 && head_of_k_maximum3!=nullptr ){
            set_same_values(head_of_k_maximum,head_of_k_maximum3);
            head_of_k_maximum->right=new Node();
            head_of_k_maximum->right->left=head_of_k_maximum;
            head_of_k_maximum=head_of_k_maximum->right;
            head_of_k_maximum3=head_of_k_maximum3->right;
            i++;
        }
        while(head_of_k_maximum->left!=nullptr){
            head_of_k_maximum=head_of_k_maximum->left;
        }
        /* while(head_of_k_maximum->right!=nullptr){
            cout<<head_of_k_maximum->book_code<<endl;
            head_of_k_maximum=head_of_k_maximum->right;
        } */
        std::cout<<".............."<<question<<endl;
        
        query_llm(filename, head_of_k_maximum, 7, "hello", gptquerry);
    }
    return;
}

void QNA_tool::query_llm(string filename, Node *root, int k, string API_KEY,string question) {
    // first write the k paragraphs into different files

    Node *traverse = root;
    int num_paragraph = 0;

    while (num_paragraph < k) {
        assert(traverse != nullptr);
        string p_file = "paragraph_";
        p_file += to_string(num_paragraph);
        p_file += ".txt";
        // delete the file if it exists
        remove(p_file.c_str());
        ofstream outfile(p_file);
        string paragraph = get_paragraph(traverse->book_code, traverse->page,
                                         traverse->paragraph);
        assert(paragraph != "$I$N$V$A$L$I$D$");
        outfile << paragraph;
        //cout << paragraph << "\n\n";
        outfile.close();
        traverse = traverse->right;
        num_paragraph++;
    }
    cout<<"\x1B[37m"<<"___________________________________________"<<"\x1B[0m"<<endl;
    
    std::cout << "\x1B[34m" << "CHAT_GPT Querry\n" << endl; 
  
    // write the query to query.txt
    ofstream outfile("query.txt");

    outfile << "These are the excerpts from Mahatma Gandhi's books.\nOn the "
               "basis of this,";
    outfile << question;
    // You can add anything here - show all your creativity and skills of using
    // ChatGPT
    outfile<<".\n There are certain rules you should follow while answering the question above: 1.Donot give vague answers like \"I donot know\",\"it is not clear\",\"the above context is not relevant\",\"it is not mentioned \" etc, 2.Answer to your best even if the context is not sufficient. 3.Try to give subjective answers 4.answer should be at least 100 words. 5.give definitive answers";
    
    outfile.close();
    cout << "These are the excerpts from Mahatma Gandhi's books.\nOn the "
               "basis of this,";
    cout<< question;
    // You can add anything here - show all your creativity and skills of using
    // ChatGPT
    cout<<".\n There are certain rules you should follow while answering the question above: 1.Donot give vague answers like \"I donot know\",\"it is not clear\",\"the above context is not relevant\",\"it is not mentioned \" etc, 2.Answer to your best even if the context is not sufficient. 3.Try to give subjective answers 4.answer should be at least 100 words. 5.give definitive answers"<<"\x1B[0m";
    // you do not need to necessarily provide k paragraphs - can configure
    // yourself

    // python3 <filename> API_KEY num_paragraphs query.txt
    cout << "\n";
    string command = "python3 ";
    command += filename;
    command += " ";
    command += API_KEY;
    command += " ";
    command += to_string(k);
    command += " ";
    command += "query.txt";
    
    std::cout << "\x1B[33m" << "CHAT_GPT Response\n" << "\x1B[0m"<<endl;
    cout<<endl;
    system(command.c_str());
    cout<<"\x1B[37m"<<"___________________________________________"<<"\x1B[0m"<<endl;
    
    return;
}

void QNA_tool::insert_sentence(int book_code, int page, int paragraph, int sentence_no, string sentence) {
    // checking if trie already exists or not

    while (book_code > StemmedParaTrie.size()) {
        vector<vector<Trie *>> n;
        ParaTrie.push_back(n);

        StemmedParaTrie.push_back(n);  // stemmed trie
    }
    while (page + 1 > StemmedParaTrie[book_code - 1].size()) {
        vector<Trie *> n;
        ParaTrie[book_code - 1].push_back(n);

        StemmedParaTrie[book_code - 1].push_back(n);  // stemmed trie
    }
    while (paragraph + 1 > StemmedParaTrie[book_code - 1][page].size()) {
        Trie *n = new Trie();
        ParaTrie[book_code - 1][page].push_back(n);

        StemmedParaTrie[book_code - 1][page].push_back(n);  // stemmed trie
    }

    // inserting in both corpus and para tries

    string word = "";

    for (int i = 0; i < sentence.size(); i++) {
        if (is_separator(sentence[i])) {
            if (word != "") {
                word = to_lower(word);  // word to lower case
                corpus->increase_freq(word);
                ParaTrie[book_code - 1][page][paragraph]->increase_freq(word);

                StemmedParaTrie[book_code - 1][page][paragraph]->increase_freq(stemmedword(word));  // stemmed trie
                StemmedCorpus->increase_freq(stemmedword(word));                                    // stemmed  trie
            }
            word = "";
        } else {
            word += sentence[i];
        }
    }

    if (word != "") {
        word = to_lower(word);
        corpus->increase_freq(word);
        ParaTrie[book_code - 1][page][paragraph]->increase_freq(word);

        StemmedParaTrie[book_code - 1][page][paragraph]->increase_freq(stemmedword(word));                           // stemmed trie
        StemmedCorpus->increase_freq(stemmedword(word));                                                             // stemmed trie
    }
}

std::string QNA_tool::get_paragraph(int book_code, int page, int paragraph) {
    cout << "Book_code: " << book_code << " Page: " << page
         << " Paragraph: " << paragraph << endl;

    std::string filename = "mahatma-gandhi-collected-works-volume-";
    filename += to_string(book_code);
    filename += ".txt";

    std::ifstream inputFile(filename);

    std::string tuple;
    std::string sentence;

    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open the input file " << filename << "."
                  << std::endl;
        exit(1);
    }

    std::string res = "";

    while (std::getline(inputFile, tuple, ')') &&
           std::getline(inputFile, sentence)) {
        // Get a line in the sentence
        tuple += ')';

        int metadata[5];
        std::istringstream iss(tuple);

        // Temporary variables for parsing
        std::string token;

        // Ignore the first character (the opening parenthesis)
        iss.ignore(1);

        // Parse and convert the elements to integers
        int idx = 0;
        while (std::getline(iss, token, ',')) {
            // Trim leading and trailing white spaces
            size_t start = token.find_first_not_of(" ");
            size_t end = token.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                token = token.substr(start, end - start + 1);
            }

            // Check if the element is a number or a string
            if (token[0] == '\'') {
                // Remove the single quotes and convert to integer
                int num = std::stoi(token.substr(1, token.length() - 2));
                metadata[idx] = num;
            } else {
                // Convert the element to integer
                int num = std::stoi(token);
                metadata[idx] = num;
            }
            idx++;
        }

        if ((metadata[0] == book_code) && (metadata[1] == page) &&
            (metadata[2] == paragraph)) {
            res += sentence;
        }
    }

    inputFile.close();
    return res;
}
