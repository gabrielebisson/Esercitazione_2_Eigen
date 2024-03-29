#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

int main()
{
    Matrix<double,2,2> A1;
    A1<< 5.547001962252291e-01,-3.770900990025203e-02, 8.320502943378437e-01, -9.992887623566787e-01;
    Matrix<double,2,1> b1;
    b1<< -5.169911863249772e-01, 1.672384680188350e-01;
    Matrix<double,2,2> A2;
    A2<< 5.547001962252291e-01, -5.540607316466765e-01, 8.320502943378437e-01, -8.324762492991313e-01;
    Matrix<double,2,1> b2;
    b2<< -6.394645785530173e-04, 4.259549612877223e-04;
    Matrix<double,2,2> A3;
    A3<< 5.547001962252291e-01, -5.547001955851905e-01, 8.320502943378437e-01,-8.320502947645361e-01;
    Matrix<double,2,1> b3;
    b3<< -6.400391328043042e-10, 4.266924591433963e-10;
    const unsigned int n=3;
    MatrixXd matrici[n]={A1,A2,A3}; //array con le matrici dei coefficienti
    MatrixXd vettori[n]={b1,b2,b3}; //array con i vettori dei termini noti (servono dopo nel ciclo)
    Matrix<double,2,1> soluz;
    soluz<< -1., -1.; //soluzione esatta

    for(unsigned int i=0;i<n;i++)
    {
        if (matrici[i].row(1).size()!=vettori[i].size()) //controlla che le dimensioni della matrice A e il vettore b siano coerenti nel sistema Ax=b
        {
            cerr<<"Non si puo' risolvere questo sistema lineare, perche' c'e' incompatibilita' tra il numero di colonne della matrice dei coefficienti ("<<matrici[i].row(1).size()<<") e la dimensione del vettore dei termini noti ("<<vettori[i].size()<<")"<<endl;
        }
        else
        {
            MatrixXd x_palu=matrici[i].lu().solve(vettori[i]);
            MatrixXd x_qr=matrici[i].householderQr().solve(vettori[i]);
            cout<<"Il sistema che ha come matrice dei coefficienti A:"<<endl<<matrici[i]<<endl;
            cout<<"e vettore dei termini noti b:"<<endl<<vettori[i]<<endl;
            cout<<"ha soluzione x_palu:"<<endl<<x_palu<<endl<<"ottenuta con la fattorizzazione della matrice A: PA=LU"<<endl;
            cout<<"con errore relativo rispetto alla soluzione esatta:"<<endl<<soluz<<endl<<"che e' pari a:"<<endl<<(x_palu-soluz).norm() / soluz.norm()<< endl;
            cout<<"Tale sistema ha inoltre soluzione x_qr:"<<endl<<x_qr<<endl<<"ottenuta con la fattorizzazione della matrice A: A=QR"<<endl;
            cout<<"con errore relativo rispetto alla soluzione esatta:"<<endl<<soluz<<endl<<"che e' pari a :"<<endl<<(x_qr-soluz).norm() / soluz.norm()<< endl;
        }
        cout<<"---------------------------------------------------"<<endl;
    }
    return 0;
}
