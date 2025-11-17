#include <petscksp.h>
#include <iostream>

/*
Questo codice usa le matrici estrapolate da Matlab e fa il calcolo. Controlla che il calcolo sia corretto
In fondo ci sono le varie norme di quanto i vari valori si discostano dalla soluzione di matlab e la exact solution
*/

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);
    /*
      =======================================================================
                      START OF MATRIX AND VECTOR BUILDING
      =======================================================================
    */
    Mat A;
    Vec b, uhm, uex; // uhm soluzione di matlab. uex exact solution

    // Open and load the matrix
    PetscViewer viewerA;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "FilesLap/matrix_A.dat", FILE_MODE_READ, &viewerA);

    MatCreate(PETSC_COMM_WORLD, &A);       // <--- CREA IL CONTENITORE PRIMA
    MatLoad(A, viewerA);                   // <--- POI CARICA I DATI
    PetscViewerDestroy(&viewerA);

    // Open and load the vector
    PetscViewer viewerB;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "FilesLap/vector_b.dat", FILE_MODE_READ, &viewerB);

    VecCreate(PETSC_COMM_WORLD, &b);       // <--- CREA IL CONTENITORE PRIMA
    VecLoad(b, viewerB);
    PetscViewerDestroy(&viewerB);

    PetscViewer viewerC;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "FilesLap/vector_uhm.dat", FILE_MODE_READ, &viewerC);

    VecCreate(PETSC_COMM_WORLD, &uhm);       // <--- CREA IL CONTENITORE PRIMA
    VecLoad(uhm, viewerC);
    PetscViewerDestroy(&viewerC);

    PetscViewer viewerD;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, "FilesLap/vector_uex.dat", FILE_MODE_READ, &viewerD);

    VecCreate(PETSC_COMM_WORLD, &uex);       // <--- CREA IL CONTENITORE PRIMA
    VecLoad(uex, viewerD);
    PetscViewerDestroy(&viewerD);


    /*
      =======================================================================
                        EDN OF MATRIX AND VECTOR BUILDING
      ======================================================================= 
    */

    Vec uh, err_1, err_2, err_3; //err_1 = uex - uh; err_2 = uh - uhm; err_3 = uex - uhm
    KSP ksp;

    VecDuplicate(b, &uh);
    VecDuplicate(b, &err_1);
    VecDuplicate(b, &err_2);
    VecDuplicate(b, &err_3);
    VecSet(uh, 0.0);


    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetTolerances(ksp, -1.0, 1e-10, PETSC_DEFAULT, PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    KSPSetType(ksp, KSPCG); //gradiente coniugato come solver

    PC pc;
    KSPGetPC(ksp, &pc); // Dice "lega il precond di ksp a pc"
    //PCSetType(pc, PCGAMG); // Setta il precondizionatore a default AMG
    //KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); // dettagli del Krilov method
    PCSetType(pc, PCHYPRE);
    PCHYPRESetType(pc, "boomeramg");  // Precond: AMG Hypre
    KSPSolve(ksp, b, uh);

    PetscReal rnorm;
    KSPGetResidualNorm(ksp, &rnorm);
    PetscPrintf(PETSC_COMM_WORLD, "Final residual norm = %.6e\n", (double)rnorm);



    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    PetscPrintf(PETSC_COMM_WORLD, "Converged in %d iterations\n", (int)its);

    VecWAXPY(err_1, -1.0, uex, uh);  // Serve per fare le operazioni tra vettori del tipo uex - 1.0 * uh
    VecWAXPY(err_2, -1.0, uex, uhm);
    VecWAXPY(err_3, -1.0, uh, uhm);

    PetscReal e1, e2, e3;
    VecNorm(err_1, NORM_2, &e1);
    VecNorm(err_2, NORM_2, &e2);
    VecNorm(err_3, NORM_2, &e3);

    std::cout << "|| uex - uh||_2 = " << e1 << std::endl;
    std::cout << "|| uex - uhm||_2 = " << e2 << std::endl;
    std::cout << "|| uhm - uh||_2 = " << e3 << std::endl;
    // The solution conververges in one step because LU factorization is not an iterative method, but a direct one
    // so we are directly solving the system, and not finding a possible one, and at each iteration we change it
    VecDestroy(&uh);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);


    PetscFinalize();
    return 0;
}
