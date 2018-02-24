//
//  tp2.cpp
//  Exemple de convolution d'image avec lodepng
//
//  Créé par Julien-Charles Lévesque
//  Copyright 2015 Université Laval. Tous droits réservés.
//


#include "lodepng.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <omp.h>

#include "Chrono.hpp"
#include "PACC/Tokenizer.hpp"

using namespace std;

//Aide pour le programme
void usage(char* inName) {
    cout << endl << "Utilisation> " << inName << " fichier_image fichier_noyau [fichier_sortie=output.png]" << endl;
    exit(1);
}

//Décoder à partir du disque dans un vecteur de pixels bruts en un seul appel de fonction
void decode(const char* inFilename,  vector<unsigned char>& outImage, unsigned int& outWidth, unsigned int& outHeight)
{
    //Décoder
    unsigned int lError = lodepng::decode(outImage, outWidth, outHeight, inFilename);

    //Montrer l'erreur s'il y en a une.
    if(lError)
        cout << "Erreur de décodage " << lError << ": " << lodepng_error_text(lError) << endl;

    //Les pixels sont maintenant dans le vecteur outImage, 4 octets par pixel, organisés RGBARGBA...
}

//Encoder à partir de pixels bruts sur le disque en un seul appel de fonction
//L'argument inImage contient inWidth * inHeight pixels RGBA ou inWidth * inHeight * 4 octets
void encode(const char* inFilename, vector<unsigned char>& inImage, unsigned int inWidth, unsigned int inHeight)
{
    //Encoder l'image
    unsigned lError = lodepng::encode(inFilename, inImage, inWidth, inHeight);

    //Montrer l'erreur s'il y en a une.
    if(lError)
        cout << "Erreur d'encodage " << lError << ": "<< lodepng_error_text(lError) << endl;
}

void executerOmp1(int inArgc, char *inArgv[]){

    // Déclaration des variables
    string lOutFilename;
    //Variables contenant des indices
    int fy, fx;
    //Variables temporaires pour les canaux de l'image
    double lR, lG, lB;

    Chrono lChrono(true);

#pragma omp parallel num_threads(4)
    {
        // cout << "Bonjour du thread " << omp_get_thread_num() << " sur " << omp_get_num_threads() << endl;
        if(inArgc < 3 or inArgc > 4) {
            usage(inArgv[0]);
        }
        string lFilename = inArgv[1];

        if (inArgc == 4)
            lOutFilename = inArgv[3];
        else
            lOutFilename = "output.png";

        // Lire le noyau.
        ifstream lConfig;
        lConfig.open(inArgv[2]);
        if (!lConfig.is_open()) {
            cerr << "Le fichier noyau fourni (" << inArgv[2] << ") est invalide." << endl;
            exit(1);
        }

        PACC::Tokenizer lTok(lConfig);
        lTok.setDelimiters(" \n","");

        string lToken;
        lTok.getNextToken(lToken);

        int lK = atoi(lToken.c_str());
        int lHalfK = lK/2;

        cout << "Taille du noyau: " <<  lK << endl;

        //Lecture du filtre
        double* lFilter = new double[lK*lK];

        for (int i = 0; i < lK; i++) {
            for (int j = 0; j < lK; j++) {
                lTok.getNextToken(lToken);
                lFilter[i*lK+j] = atof(lToken.c_str());
            }
        }

        //Lecture de l'image
        //Variables à remplir
        unsigned int lWidth, lHeight;
        vector<unsigned char> lImage; //Les pixels bruts
        //Appeler lodepng
        decode(lFilename.c_str(), lImage, lWidth, lHeight);

        // #pragma omp for schedule(static, 4)
        for(int x = lHalfK; x < (int)lWidth - lHalfK; x++)
        {
            for (int y = lHalfK; y < (int)lHeight - lHalfK; y++)
            {
                lR = 0.;
                lG = 0.;
                lB = 0.;
                for (int j = -lHalfK; j <= lHalfK; j++) {
                    fy = j + lHalfK;
                    for (int i = -lHalfK; i <= lHalfK; i++) {
                        fx = i + lHalfK;
                        //R[x + i, y + j] = Im[x + i, y + j].R * Filter[i, j]
                        lR += double(lImage[(y + j)*lWidth*4 + (x + i)*4]) * lFilter[fx + fy*lK];
                        lG += double(lImage[(y + j)*lWidth*4 + (x + i)*4 + 1]) * lFilter[fx + fy*lK];
                        lB += double(lImage[(y + j)*lWidth*4 + (x + i)*4 + 2]) * lFilter[fx + fy*lK];
                    }
                }
                //Placer le résultat dans l'image.
                lImage[y*lWidth*4 + x*4] = (unsigned char)lR;
                lImage[y*lWidth*4 + x*4 + 1] = (unsigned char)lG;
                lImage[y*lWidth*4 + x*4 + 2] = (unsigned char)lB;
            }
        }

        //Sauvegarde de l'image dans un fichier sortie
        encode(lOutFilename.c_str(),  lImage, lWidth, lHeight);

        delete lFilter;
    };

    cout << "L'image a été filtrée et enregistrée dans " << lOutFilename << " avec succès!" << endl;

    cout << "Temp d'exécution omp 1:  " << lChrono.get() << endl;

    lChrono.pause();

}

void executerSequentiel(int inArgc, char *inArgv[]){

    Chrono lChrono(true);

    if(inArgc < 3 or inArgc > 4) usage(inArgv[0]);
    string lFilename = inArgv[1];
    string lOutFilename;
    if (inArgc == 4)
        lOutFilename = inArgv[3];
    else
        lOutFilename = "output.png";

    // Lire le noyau.
    ifstream lConfig;
    lConfig.open(inArgv[2]);
    if (!lConfig.is_open()) {
        cerr << "Le fichier noyau fourni (" << inArgv[2] << ") est invalide." << endl;
        exit(1);
    }

    PACC::Tokenizer lTok(lConfig);
    lTok.setDelimiters(" \n","");

    string lToken;
    lTok.getNextToken(lToken);

    int lK = atoi(lToken.c_str());
    int lHalfK = lK/2;

    cout << "Taille du noyau: " <<  lK << endl;

    //Lecture du filtre
    double* lFilter = new double[lK*lK];

    for (int i = 0; i < lK; i++) {
        for (int j = 0; j < lK; j++) {
            lTok.getNextToken(lToken);
            lFilter[i*lK+j] = atof(lToken.c_str());
        }
    }

    //Lecture de l'image
    //Variables à remplir
    unsigned int lWidth, lHeight;
    vector<unsigned char> lImage; //Les pixels bruts
    //Appeler lodepng
    decode(lFilename.c_str(), lImage, lWidth, lHeight);

    //Variables contenant des indices
    int fy, fx;
    //Variables temporaires pour les canaux de l'image
    double lR, lG, lB;
    for(int x = lHalfK; x < (int)lWidth - lHalfK; x++)
    {
        for (int y = lHalfK; y < (int)lHeight - lHalfK; y++)
        {
            lR = 0.;
            lG = 0.;
            lB = 0.;
            for (int j = -lHalfK; j <= lHalfK; j++) {
                fy = j + lHalfK;
                for (int i = -lHalfK; i <= lHalfK; i++) {
                    fx = i + lHalfK;
                    //R[x + i, y + j] = Im[x + i, y + j].R * Filter[i, j]
                    lR += double(lImage[(y + j)*lWidth*4 + (x + i)*4]) * lFilter[fx + fy*lK];
                    lG += double(lImage[(y + j)*lWidth*4 + (x + i)*4 + 1]) * lFilter[fx + fy*lK];
                    lB += double(lImage[(y + j)*lWidth*4 + (x + i)*4 + 2]) * lFilter[fx + fy*lK];
                }
            }
            //Placer le résultat dans l'image.
            lImage[y*lWidth*4 + x*4] = (unsigned char)lR;
            lImage[y*lWidth*4 + x*4 + 1] = (unsigned char)lG;
            lImage[y*lWidth*4 + x*4 + 2] = (unsigned char)lB;
        }
    }

    //Sauvegarde de l'image dans un fichier sortie
    encode(lOutFilename.c_str(),  lImage, lWidth, lHeight);

    lChrono.pause();

    cout << "L'image a été filtrée et enregistrée dans " << lOutFilename << " avec succès!" << endl;

    cout << "Temp d'exécution : " << lChrono.get() << endl;

    delete lFilter;
}

// Code Diego
void executerParallele(int inArgc, char *inArgv[]){

    if(inArgc < 3 or inArgc > 4) {
        usage(inArgv[0]);
    }
    string lFilename = inArgv[1];
    string lOutFilename;
    if (inArgc == 4)
        lOutFilename = inArgv[3];
    else
        lOutFilename = "output.png";

    // Lire le noyau.
    ifstream lConfig;

    lConfig.open(inArgv[2]);
    if (!lConfig.is_open()) {
        cerr << "Le fichier noyau fourni (" << inArgv[2] << ") est invalide." << endl;
        exit(1);
    }

    PACC::Tokenizer lTok(lConfig);
    lTok.setDelimiters(" \n","");

    string lToken;
    lTok.getNextToken(lToken);

    int lK = atoi(lToken.c_str());
    int lHalfK = lK/2;

    Chrono lChrono(true);

    //Lecture du filtre
    double* lFilter;
    lFilter = new double[lK * lK];

    //Lecture de l'image
    //Variables à remplir
    unsigned int lWidth, lHeight;
    vector<unsigned char> lImage; //Les pixels bruts

    // Doesnt matter what we try here, it breaks the image sometimes
    // #pragma omp parallel for // private(lToken) // schedule(auto) // collapse(2)
#pragma omp for ordered
    for (int i = 0; i < lK; i++) {
        for (int j = 0; j < lK; j++) {
            //#pragma omp atomic read
            lTok.getNextToken(lToken);
            lFilter[i * lK + j] = atof(lToken.c_str());
        }
    }

    //Appeler lodepng
    decode(lFilename.c_str(), lImage, lWidth, lHeight);

    //Variables contenant des indices
    int fy, fx;
    //Variables temporaires pour les canaux de l'image
    double lR, lG, lB;

    int maxWidth = (int)lWidth - lHalfK;
    int maxHeight = (int)lHeight - lHalfK;

    //omp_set_num_threads(maxWidth);

    //#pragma omp parallel for schedule(dynamic) shared(lImage, lWidth, lHeight) private(lR, lG, lB, fy, fx)

    //#pragma omp parallel for schedule(dynamic) shared(lImage, lWidth, lHeight) private(lR, lG, lB, fy, fx)

    // collapse(2) breaks sometimes the image
#pragma omp parallel for schedule(static, 4) shared(lImage, lWidth, lHeight) private(lR, lG, lB, fy, fx)
    for(int x = lHalfK; x < maxWidth; x++)
    {
        for (int y = lHalfK; y < maxHeight; y++)
        {
//#pragma omp atomic write
            lR = 0.;
//#pragma omp atomic write
            lG = 0.;
//#pragma omp atomic write
            lB = 0.;

            //This does not work at all - 40 sec +/-
            //#pragma omp parallel for schedule(dynamic) shared(lImage, lWidth, lHeight)
            // #pragma omp for ordered
            for (int j = -lHalfK; j <= lHalfK; j++) {
                fy = j + lHalfK;
                for (int i = -lHalfK; i <= lHalfK; i++) {
                    fx = i + lHalfK;

                    //#pragma omp atomic update
                    lR += double(lImage[(y + j)*lWidth*4 + (x + i)*4]) * lFilter[fx + fy*lK];
                    //#pragma omp atomic update
                    lG += double(lImage[(y + j)*lWidth*4 + (x + i)*4 + 1]) * lFilter[fx + fy*lK];
                    //#pragma omp atomic update
                    lB += double(lImage[(y + j)*lWidth*4 + (x + i)*4 + 2]) * lFilter[fx + fy*lK];

                }
            }
            //Placer le résultat dans l'image.
            lImage[y*lWidth*4 + x*4] = (unsigned char)lR;
            lImage[y*lWidth*4 + x*4 + 1] = (unsigned char)lG;
            lImage[y*lWidth*4 + x*4 + 2] = (unsigned char)lB;
        }
    }

    //Sauvegarde de l'image dans un fichier sortie
    encode(lOutFilename.c_str(),  lImage, lWidth, lHeight);

    delete lFilter;

    lChrono.pause();

    cout << "L'image a été filtrée et enregistrée dans " << lOutFilename << " avec succès!" << endl;

    cout << "Temps d'execution parallele  = \033[1;31m" << lChrono.get() << " sec\033[0m" << endl;

}

int main(int inArgc, char *inArgv[])
{
    cout << "Execution séquentielle " << endl;
    // executerSequentiel(inArgc, inArgv);
    cout << "Execution OMP 1" << endl;
    //executerOmp1(inArgc, inArgv);
    cout << "Execution OMP 2" << endl;
    executerParallele(inArgc, inArgv);
    return 0;
}