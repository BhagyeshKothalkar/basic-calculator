#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
using namespace std;
const double e = 2.718281828459045;
double ln2 = 0.693147180559945;
const double prdctcvalues = 0.607252935008881;
const double pi = 3.14159265358979;
    double arctanvalues[51] = {
        0.785398163397448,  // arctan(1/2^0)
        0.463647609000806,  // arctan(1/2^1)
        0.244978663126864,  // arctan(1/2^2)
        0.124354994546761,  // arctan(1/2^3)
        0.062418809995957,  // arctan(1/2^4)
        0.031239833430268,  // arctan(1/2^5)
        0.015623728620477,  // arctan(1/2^6)
        0.007812341060101,  // arctan(1/2^7)
        0.003906230131967,  // arctan(1/2^8)
        0.001953122516479,  // arctan(1/2^9)
        0.000976562189559,  // arctan(1/2^10)
        0.000488281211195,  // arctan(1/2^11)
        0.000244140620149,  // arctan(1/2^12)
        0.000122070311893,  // arctan(1/2^13)
        0.000061035156174,  // arctan(1/2^14)
        0.000030517578115,  // arctan(1/2^15)
        0.000015258789061,  // arctan(1/2^16)
        0.000007629394531,  // arctan(1/2^17)
        0.000003814697266,  // arctan(1/2^18)
        0.000001907348633,  // arctan(1/2^19)
        0.000000953674316,  // arctan(1/2^20)
        0.000000476837158,  // arctan(1/2^21)
        0.000000238418579,  // arctan(1/2^22)
        0.000000119209290,  // arctan(1/2^23)
        0.000000059604645,  // arctan(1/2^24)
        0.000000029802322,  // arctan(1/2^25)
        0.000000014901161,  // arctan(1/2^26)
        0.000000007450581,  // arctan(1/2^27)
        0.000000003725290,  // arctan(1/2^28)
        0.000000001862645,  // arctan(1/2^29)
        0.000000000931323,  // arctan(1/2^30)
        0.000000000465661,  // arctan(1/2^31)
        0.000000000232831,  // arctan(1/2^32)
        0.000000000116415,  // arctan(1/2^33)
        0.000000000058208,  // arctan(1/2^34)
        0.000000000029104,  // arctan(1/2^35)
        0.000000000014552,  // arctan(1/2^36)
        0.000000000007276,  // arctan(1/2^37)
        0.000000000003638,  // arctan(1/2^38)
        0.000000000001819,  // arctan(1/2^39)
        0.000000000000909,  // arctan(1/2^40)
        0.000000000000454,  // arctan(1/2^41)
        0.000000000000227,  // arctan(1/2^42)
        0.000000000000113,  // arctan(1/2^43)
        0.000000000000057,  // arctan(1/2^44)
        0.000000000000028,  // arctan(1/2^45)
        0.000000000000014,  // arctan(1/2^46)
        0.000000000000007,  // arctan(1/2^47)
        0.000000000000003,  // arctan(1/2^48)
        0.000000000000002,  // arctan(1/2^49)
        0.000000000000001   // arctan(1/2^50)
    };

    

void add()
{
    double x,y;
    cout<<"\nEnter the two numbers: ";
    cin>>x>>y;
    cout<<"sum of "<<setprecision(15)<<x<<" and "<<setprecision(15)<<y<<" is: "<<setprecision(15)<<x+y;
}
void subtract()
{
    double x,y;
    cout<<"\nEnter the two numbers: ";
    cin>>x>>y;
    cout<<"difference between "<<setprecision(15)<<x<<"-"<<setprecision(15)<<y<<" is: "<<setprecision(15)<<x+y;
}
void multiply()
{
    double x,y;
    cout<<"\nEnter the two numbers: ";
    cin>>x>>y;
   cout<<"product of "<<setprecision(15)<<x<<" and "<<setprecision(15)<<y<<" is: "<<setprecision(15)<<x+y;
}
void divide()
{
    double x,y;
    cout<<"\nEnter the two numbers: ";
    cin>>x>>y;
    cout<<"value of "<<setprecision(15)<<x<<"/"<<setprecision(15)<<y<<" is: "<<setprecision(15)<<x+y;
}
double raiseto(double base, int power)//raises a double base to an integer power
{
    double answer = 1;
    for(int i = 0; i< power; i++){
        answer*= base;
    }
    return(answer);
}
double sqrt (double a){
    double s = 1;int i = 0;
    while (i < 50){
        s = (s + a/s)/2;
        i++;
    }
    return s;
}

void solvequad()//asks for coeffecients ofquadrtic equation and prints its roots 
{
    double a,b,c;
    cout<<"\nEnter the coeffecient of (x^2): ";
    cin>>a;cout<<a;
    cout<<"\nEnter the coeffecient of (x): ";
    cin>>b;cout<<b;
    cout<<"\nEnter the constant term: ";
    cin>>c;cout<<c;
    double discriminant = b*b - 4*a*c;
    cout<<"discrminant :"<<discriminant<<endl;
    if (discriminant > 0)//if discriminant is positive then roots are real and distinct. thus we calculate using quadratic formula
    {
        cout<<"\n roots are: "<<setprecision(15)<<(-b + sqrt(discriminant))/(2*a)<<" and "<<setprecision(15)<<(-b - sqrt(discriminant))/(2*a); 
    }
    else if (discriminant == 0)//if discriminant is zero then roots are real and equal
    {
        cout<<"\n roots are equal, and are equal to: "<<-b/(2*a); 
    }
    else//if discriminant is negative then roots are non real complex. currently these are not supported.
    {cout<<"the roots are non real complex";}
}
double ln(double z = 0, bool mode = 1){
    if(mode){
    cout<<"enter the number: ";
    double z; cin>>z;}
    double divisions = 0;
    while (z > 1)
    {
        z /= 2;
        divisions++;
    }
    while (z < 0.5)
    {
        z*= 2; 
        divisions --;
    }
    //to calculate log by taylor series, we use ln (1 - x) thus, the ln(z) = ln(1-(1-z))
    //and this 1-z is the one whose series eill be formed. 
    z = 1-z;
    double term = z;
    double sum = z;
    for(int i = 2; i < 20; i++){//here, i represents the i th term. 
    //we have initialized the first term and now we go for the other terms.
        term *= z*(i-1)/(i);
        sum += term;
    }
    if(mode){cout<<setprecision(15)<<divisions*ln2-sum<<endl;}
    else{return divisions*ln2-sum;}
}
double exp(double x){
    int y = x;
    double out = raiseto(e,y);
    double z = x-y;
    double term = 1,sum = 1;
    for(int i = 1; i < 50; i++){
        term *= (z/i);
        sum += term;
    }
    return out*sum;
}
void anyexp(){
    cout<<"\n enter the base: ";
    double base; cin >> base;
    cout<<"\nenter exponent: ";
    double exponent; cin>>exponent;
    cout<<exp(exponent*ln(base,false));
}
double determinant(vector<vector<double>> matrix)
{
    int dimension = matrix.size();
    double out = 0;
    if (dimension == 2){
        return (matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1]);
        }
    else{
    for(int i = 0; i < dimension; i++)
    {
        vector<vector<double>> mat = matrix;
        double element = mat[0][i];
        mat.erase(mat.begin());
        for(int j = 0; j < dimension - 1; j++){
            mat[j].erase(mat[j].begin() + i);
        }
        if(i % 2 == 0) {out += (element*determinant(mat));}
        else {out -= (element*determinant(mat));}
    }
    return out;
    }
    
}

void solveeq(double detcoeffmatrix, vector<double> & detansmatrix){
    if(detcoeffmatrix !=0 ){
        for(int i = 0; i < detansmatrix.size(); i++)
        {
           detansmatrix[i] /= detcoeffmatrix; 
        }
    }
}
void entervalues() {
    string n = ".0123456789";
    string s = "abcdefghijklnmopqrstuvwxyz";
    cout<<"\nenter number of variables(also equal to number of equations): ";
    int varno; cin>>varno;
    vector<vector<double>> coeffmatrix(varno, vector<double>(varno));
    vector<char>varmatrix(varno);
    vector<double> ansmatrix(varno),detansmatrix(varno);
    int i = 0;//here, i and j are to represent the columnnumber and rownumber 
                    //of the coeffecient matrix we eventually fill

    for(int j = 0; j < varno; j++) {
        cout << "\nenter the equation  "<<j+1;
        cout<<"\n  or enter 'listout' for the conventions to be followed to avoid any error"<<endl;
        
        string equation;
        cin >> equation;
        if(equation == "listout"){
            cout<<"*the equations should consist of all variables on one side of equality sign and the constant terms on the other";
            cout<<"\n*each equaton should consist of all variables involved, with some numeric coeffecients";
            cout<<"\n*(that means that any variable, (say a) not appearinging an equation has to be shown as 0a)";
            cout<<"\n*also, the variable (say, b) occuring with coeffecient 1 has to be written along with the coeffecient(1b cannot be written as b)";
            cout<<"\n*the coeffecients should strictly be real numbers(positive or negative)and not expressions.";
            cout<<"\n*the signs + and - cannot be next to each other in the expression";
            cout<<"\n*the order of variables used in equation cannot be changed for example, if variables occur in the equation 1 as 'a' 'b' 'c'";
            cout<<"\n*then they should occur in same order in all equations that follow";
            cout<<"\n*e.g. 23.5a - 69.65b = -23.65 is a valid equation";
            --j;continue;
        }
        
        string tovec;
        for (int k = 0; k < equation.length(); k++)//k is the index of the character in the string over which we iterate
        {
            if (n.find(equation[k]) != string::npos) {
                tovec += equation[k];
            } 
            else if (s.find(equation[k]) != string::npos) {
                if (i < coeffmatrix.size() && i < coeffmatrix[j].size()) {
                    coeffmatrix[j][i] = stod(tovec);
                    tovec = "";
                    varmatrix[i] = equation[k];
                    i++; 
                }
                
            } 
            else if (equation[k] == '-') {
                tovec = "-";
            } 
            else if (equation[k] == '+') {
                tovec = "";
            } 
            else if (equation[k] == '=') {
                tovec = "-";
            }
        }
        if (j < varno) {
            ansmatrix[j] = -stod(tovec);
        }
        i = 0;
    }
    double detcoeffmatrix = determinant(coeffmatrix);
    vector<vector<double>> coeffmatrix1(varno, vector<double>(varno));
    for(int i = 0; i < varno; i++)
    {
        coeffmatrix1 = coeffmatrix;
        for(int j = 0; j < varno; j++)
        {
            coeffmatrix1[j][i] = ansmatrix[j];
        }
        detansmatrix[i] = determinant(coeffmatrix1);
    }
    //coutarray(coeffmatrix);
    for(int i = 0; i < varno; i++){
        detansmatrix[i] /= detcoeffmatrix; 
        cout<<"\n the value of "<<varmatrix[i]<<" is: "<<setprecision(15)<<detansmatrix[i];
    }
}
void CORDICR(double& x, double& y,double& z, double refarr[],int size)
/*runs CORDIC algorithm in ROTATION MODE. reference: https://people.clas.ufl.edu/bruceedwards/files/paper1.pdf */
{
    
    for(int i = 0; i < size; i++){
        double a = x,b = y;
        if (z < 0) {
            x = a + (b / (1 << i));
            y = b - (a / (1 << i));
            z += refarr[i];
            }
         else {
            x = a - (b / (1 << i));
            y= b + (a / (1 << i));
            z -= refarr[i];
        }
    }
    x *= prdctcvalues;
    y *= prdctcvalues;
}
void trig(string mode){
    cout<< "\nenter the angle in degrees: ";
    double angle; cin>>angle;
    double x = 1, y = 0, z = angle*(pi/180);
    if (mode == "sin")
    {   
        while(z > pi){z -= 2*pi;}
        while(z < -pi){z += 2*pi;}
        if(z > pi/2){z = pi -z;}
        else if(z < -pi/2){z = -pi -z;}
        CORDICR(x,y,z,arctanvalues,51);
        cout<<"sin of "<<setprecision(15)<<angle<<" is: "<<setprecision(15)<<y<<endl;
    }
    else if (mode == "cos")
    {
        int sign = 1;
        while(z > pi){z -= 2*pi;}
        while(z < -pi){z += 2*pi;}
        if(z > pi/2){z = pi -z;sign *= -1;}
        else if(z < -pi/2){z = -pi -z;sign *= -1;}
        CORDICR(x,y,z,arctanvalues,51);
        cout <<"cos of "<<setprecision(15)<<angle<<" is: "<<setprecision(15)<< sign*x<<endl;
    }
    else if (mode == "tan"&& x!= 0){
        while(z > pi){z -= 2*pi;}
        while(z < -pi){z += 2*pi;}
        if(z > pi/2){z = z - pi;}
        else if(z < -pi/2){z = z + pi;}
        CORDICR(x,y,z,arctanvalues,51);
        cout <<"tan of "<<setprecision(15)<<angle<<" is: "<<setprecision(15)<< y/x <<endl;
    }
    else if (mode == "cot" && y != 0){
        while(z > pi){z -= 2*pi;}
        while(z < -pi){z += 2*pi;}
        if(z > pi/2){z = pi -z;}
        else if(z < -pi/2){z = -pi -z;}
        CORDICR(x,y,z,arctanvalues,51);
        cout <<"cot of "<<setprecision(15)<<angle<<" is: "<<setprecision(15)<< x/y<<endl;
    }
    else if (mode == "sec"&& x!= 0){
        int sign = 1;
        while(z > pi){z -= 2*pi;}
        while(z < -pi){z += 2*pi;}
        if(z > pi/2){z = pi -z;sign *= -1;}
        else if(z < -pi/2){z = -pi;sign *= -1;;}
        CORDICR(x,y,z,arctanvalues,51);
        cout <<"sec of "<<setprecision(15)<<angle<<" is: "<<setprecision(15)<< sign*1/x<<endl;
    }
    else if (mode == "cosec" && y != 0){
        while(z > pi){z -= 2*pi;}
        while(z < -pi){z += 2*pi;}
        if(z > pi/2){z = pi -z;}
        else if(z < -pi/2){z = -pi -z;}
        CORDICR(x,y,z,arctanvalues,51);
        cout <<"cosec of "<<setprecision(15)<<angle<<" is: "<<setprecision(15)<< 1/y <<endl;} 
}
void CORDICV(double& x, double& y,double& z, double refarr[],int size)/*this function runs the CORDIC algorithm in VECTORING mode.
 in this mode, it could effeciently calculate inverse trigonometric and
inverse hyperboic functions. it is essentially based on rotating a unit vector having given component to to y axis,
till it becomes parallel to x axis.*/
{
    x *= prdctcvalues;//scaling down
    y *= prdctcvalues;    
    for(int i = 0; i < size; i++){
        double a = x,b = y;
        if (y > 0) {//we try to push ytowards 0 byrotation
            x = a + (b / (1 << i));
            y = b - (a / (1 << i));//if y is positive, rotate so as to decrease it.
            z += refarr[i];//as this rotation is clockwise, the angle will decrease
        }
         else{
            x = a - (b / (1 << i));
            y = b + (a / (1 << i));//if y is negative, rotate so as to increase it..
            z -= refarr[i];//as the rotation is anticlockwise, the angle will increase..
        }
    }
}
void trig1(string mode){
    cout<< "\nenter the value: ";
    double value; cin>>value;
    double x , y , z = 0;
    
    if (mode == "arcsin")
    {
        y = value; x = sqrt(1 - y*y);
        if (y >= 0){CORDICV(x,y,z,arctanvalues,29);cout<<"arcsin of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<<z;}
        else{y *= -1;CORDICV(x,y,z,arctanvalues,29);cout<<"arcsin of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< -z;}
    }
    else if (mode == "arccos")
    {
        x = value; y = sqrt(1 - x*x);
        if (x >= 0){CORDICV(x,y,z,arctanvalues,29);cout<<"arccos of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< z;}
        else{x *= -1;CORDICV(x,y,z,arctanvalues,29);cout<<"arccos of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< (pi - z);}
    }
    else if (mode == "arctan")
    {
        x = sqrt(1/(1+ value*value)); y = sqrt(1 - x*x);
        CORDICV(x,y,z,arctanvalues,29);
        if(value >= 0){cout<<"arctan of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< z;}else {cout<<"arctan of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< -z;}
    }
    else if (mode == "arccot")
    {
        y = sqrt(1/(1+ value*value)); x = sqrt(1 - y*y);
        CORDICV(x,y,z,arctanvalues,29);
        if(value >= 0){cout<<"arccot of "<<setprecision(15)<<value<<"is: " << z;} else {cout<<"arccot of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< (pi-z);}
    }
    else if (mode == "arcsec")
    {
        x = 1/value; y = sqrt(1 - x*x);
        if (x > 0){CORDICV(x,y,z,arctanvalues,29);cout<<"arccot of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< z;}
        else{x *= -1;CORDICV(x,y,z,arctanvalues,29);cout<<"arccot of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< (pi - z);}
    }
    else if (mode == "arccosec")
    {
        y = 1/value; x = sqrt(1 - y*y);
        if (y > 0){CORDICV(x,y,z,arctanvalues,29);cout<<"arccot of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< (z);}
        else{y *= -1;CORDICV(x,y,z,arctanvalues,29);cout<<"arccot of "<<setprecision(15)<<value<<"is: " <<setprecision(15)<< (z*-1);}
    }    
}

int main()
{   
    while (true)
    {
        string input;
        cout<<"\n enter the operation or enter \"listout\" for possible commands"<<"\n";
        cin>>input;
        if (input == "add"){add();}
        else if (input == "subtract"){subtract();}
        else if (input == "multiply"){multiply();}
        else if (input == "divide"){divide();}
        else if (input == "raiseto"){anyexp();}
        else if (input == "sin"){trig("sin");}
        else if (input == "cos"){trig("cos");}
        else if (input == "exp"){cout<<exp(5.6);}
        else if (input == "tan"){trig("tan");}
        else if (input == "cot"){trig("cot");}
        else if (input == "sec"){trig("sec");}
        else if (input == "cosec"){trig("cosec");}
        else if (input == "quadratic"){solvequad();}
        else if (input == "arcsin"){trig1("arcsin");}
        else if (input == "arccos"){trig1("arccos");}
        else if (input == "arctan"){trig1("arctan");}
        else if (input == "arccot"){trig1("arccot");}
        else if (input == "arcsec"){trig1("arcsec");}
        else if (input == "arccosec"){trig1("arccosec");}
        else if (input == "solvelinear"){entervalues();}
        else if (input == "over"){break;}
        else if (input == "listout"){
            cout<<"add             "<<"adds two numbers                               ";cout<<"subtract        "<<"subtract two numbers"<<endl;
            cout<<"multiply        "<<"multiply two numbers                           ";cout<<"divide          "<<"divide a number by other"<<endl;
            cout<<"raiseto         "<<"raises a base(double) to an integer power      ";cout<<"quadratic       "<<"gives root of quadratic"<<endl;
            cout<<"sin             "<<"calculates sin of given number                 ";cout<<"cos             "<<"calculates cos of given number"<<endl;
            cout<<"tan             "<<"calculates tan of given number                 ";cout<<"cot             "<<"calculates cot of given number"<<endl;
            cout<<"sec             "<<"calculates sec of given number                 ";cout<<"cosec           "<<"calculates cisec of given number"<<endl;
            cout<<"arcsin          "<<"calculates arcsin of given number              ";cout<<"arccos          "<<"calculates arccos of given number"<<endl;
            cout<<"arctan          "<<"calculates arctan of given number              ";cout<<"arccot          "<<"calculates arccot of give n number"<<endl;
            cout<<"arcsec          "<<"calculates arcsec of given number              ";cout<<"arccosec        "<<"calculates arccosec of given number"<<endl;
            cout<<"solvelinear     "<<"solves a system of linear equations            ";cout<<"over            "<<"ends the program"<<endl;
        }
        else{cout<<"your command did not match any of our operations.\nPlease try again";}
    }
    return 0;
}