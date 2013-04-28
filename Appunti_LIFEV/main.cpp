//@HEADER
/*
*******************************************************************************
VERSIONE COPIATA DA LIFEV PER PRENDERE APPUNTI!!!
*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief

    @author
    @date
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#else
#include <lifev/core/filter/ExporterVTK.hpp>
#endif

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCManage.hpp>

#include <lifev/core/mesh/RegionMesh1DStructured.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>

using namespace LifeV;


namespace
{
// Registrazione della factory ---> ??
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}

Real exactSolution ( const Real& /* t */, const Real& x, const Real& /*y*/, const Real& /* z */, const ID& /* i */ )
{
//Funzione soluzione esatta sin(pi/2x)
    return std::sin ( M_PI * 0.5 * x );
}


Real fRhs ( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */ , const ID& /* i */ )
{
//Termine noto 1.25*pi^2sin(pi/2*x)
    return 1.25 * M_PI * M_PI * std::sin ( M_PI * 0.5 * x );
}


int
main ( int argc, char* argv[] )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif

    {// a cosa serve??
        // needed to properly destroy all objects inside before mpi finalize

#ifdef HAVE_MPI
		//Puntatore a comunicatore, il comunicatore è un oggetto che contiene le varie informazioni
		//sui processori in uso.
        boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
		//Questo assert spiega che è inutile girare in parallelo se si ha un solo processore        
		ASSERT ( Comm->NumProc() < 2, "The test does not run in parallel." );
#else
        boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

		// ResgionMesh<LinearLine>, linearline è il tipo di elemento cioè un segmento
		// --> il tipo è dunque una mesh1D
        typedef RegionMesh<LinearLine> mesh_Type;
		// MatrixEpetra contenente dei reali! Chi è veramente MatrixEpetra?
        typedef MatrixEpetra<Real> matrix_Type;
        typedef VectorEpetra vector_Type;
		//Il tipo di spazio FEM dipende dalla mesh e dalla mappa.
        typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
		//Puntatore a spazio FEM
        typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

		//In questo modo solo il MASTER è verboso
        const bool verbose (Comm->MyPID() == 0);

        // Read first the data needed

        if (verbose)
        {
            std::cout << " -- Reading the data ... " << std::flush;
        }
		//Come funziona getpot??
        GetPot dataFile ( "data" );
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Build the mesh

        if (verbose)
        {
            std::cout << " -- Reading the mesh ... " << std::flush;
        }
		//Creo un oggetto di tipo MeshData, dovrebbe contenere le informazioni sul dominio più altro sulla mesh.
        MeshData meshData (dataFile, "mesh");
		//Creo un puntatore a un oggetto di tipo mesh inizializzandolo a una mesh vuota che però conosce il comunicatore
        boost::shared_ptr< mesh_Type > meshPtr ( new mesh_Type ( Comm ) );

		//Ecco qui viene effettivamente costruita la mesh
		// il primo argomento è di ritorno, il secondo BOH?
		// 20 elementi non verbosa su (0,1)
        // Set up the structured mesh
        regularMesh1D ( *meshPtr, 0,
                        dataFile ( "mesh/n", 20 ),
                        dataFile ( "mesh/verbose", false ),
                        dataFile ( "mesh/length", 1. ),
                        dataFile ( "mesh/origin", 0. ) );

        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Build the FESpaces
        if (verbose)
        {
            std::cout << " -- Building FESpaces ... " << std::flush;
        }
		//Qui viene costrutito e assegnato a un puntatore lo spazio elementi finiti
		// dipende dalla mesh dalla formula di quadratura ..??
        feSpacePtr_Type uFESpace ( new feSpace_Type ( meshPtr, feSegP1, quadRuleSeg1pt, quadRuleNode1pt, 1, Comm ) );
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }
        if (verbose)
        {
            std::cout << " ---> Dofs: " << uFESpace->dof().numTotalDof() << std::endl;
        }

        // Build the assembler and the matrices

        if (verbose)
        {
            std::cout << " -- Building assembler ... " << std::flush;
        }
		//Inizializzo l'oggetto che assemblerà la matrice
		//Gli spiego quali matrici e vettori userò e che tipo di mesh ho
        ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
        if (verbose)
        {
            std::cout << " done!" << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Setting up assembler ... " << std::flush;
        }
		//Gli assegno lo spazio elementi finiti
        adrAssembler.setFespace (uFESpace);
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Defining the matrix ... " << std::flush;
        }
		//Creo un puntatore a una matrice epetra e la inizializzo con la mappa dello spazio elemento finito
		//Quando è stata inizializzata quella mappa? fa da solo lo spazio FEM?
        boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uFESpace->map() ) );
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }

        // Perform the assembly of the matrix

        if (verbose)
        {
            std::cout << " -- Adding the diffusion ... " << std::flush;
        }
		//Aggiungo la diffusione mu=1 costante
        adrAssembler.addDiffusion (systemMatrix, 1.);
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }
        if (verbose)
        {
            std::cout << " Time needed : " << adrAssembler.diffusionAssemblyChrono().diffCumul() << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Adding the mass ... " << std::flush;
        }
		//Aggiungo la reazione sigma=pi^2
        adrAssembler.addMass (systemMatrix, M_PI * M_PI );
        if (verbose)
        {
            std::cout << " done! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Closing the matrix ... " << std::flush;
        }
		//Chiusura della matrice (??)
        systemMatrix->globalAssemble();
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Definition and assembly of the RHS
        if (verbose)
        {
            std::cout << " -- Building the RHS ... " << std::flush;
        }
		//Creo un vettore distribuito sui vari processori con un elemento ripetuto per ottimizzare le performance
        vector_Type rhs (uFESpace->map(), Repeated);
		//Ne creo un altro uguale
        vector_Type fInterpolated (uFESpace->map(), Repeated);
		//Chiamo il metodo interpolate dello spazio FEM-> come funziona??
		//Cmq il risultato lo infilo dentro fInterpolated
		uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( fRhs ), fInterpolated, 0.0 );
		//Non è chiarissimo perchè si debba chiamare addMassRhS comunque inserisce crea il vettore rhs a partire dalla forzante interpolata
		//sullo spazio FEM	
        adrAssembler.addMassRhs (rhs, fInterpolated);
        rhs.globalAssemble();//Chiudo anche il vettore (??)

        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Definition and application of the BCs
        if (verbose)
        {
            std::cout << " -- Building the BCHandler ... " << std::flush;
        }
		//Oggetto che tiene le BC
        BCHandler bchandler;

		//Creo un oggetto di tipo BCFunctionBase inizializzato con la soluzione esatta
        BCFunctionBase BCu ( static_cast<feSpace_Type::function_Type> ( exactSolution ) );
		//Aggiungo due condizioni di tipo dirichlet-->quanti parametri!!
        bchandler.addBC ("Dirichlet", Structured1DLabel::LEFT, Essential, Full, BCu, 1);
        bchandler.addBC ("Dirichlet", Structured1DLabel::RIGHT, Essential, Full, BCu, 1);

        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Updating the BCs ... " << std::flush;
        }
		//che fa?
        bchandler.bcUpdate (*uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Applying the BCs ... " << std::flush;
        }
		//Creo un vettore distribuito stavolta non repeated  (Perchè?)
        vector_Type rhsBC (rhs, Unique);
		//Applico le BC, come?? dove finisce rhs se viene sovrascritto?
        bcManage (*systemMatrix, rhsBC, *uFESpace->mesh(), uFESpace->dof(), bchandler, uFESpace->feBd(), 1.0, 0.0);
        rhs = rhsBC;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Definition of the solver
        if (verbose)
        {
            std::cout << " -- Building the solver ... " << std::flush;
        }
		//Mi sembra di ricordare che questo è da cambiare-->non val la pena impararlo
        SolverAztecOO linearSolver;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Setting up the solver ... " << std::flush;
        }
        linearSolver.setDataFromGetPot (dataFile, "solver");
        linearSolver.setupPreconditioner (dataFile, "prec");
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Setting matrix in the solver ... " << std::flush;
        }
        linearSolver.setMatrix (*systemMatrix);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        linearSolver.setCommunicator (Comm);

        // Definition of the solution
        if (verbose)
        {
            std::cout << " -- Defining the solution ... " << std::flush;
        }
        vector_Type solution (uFESpace->map(), Unique);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Solve the solution
        if (verbose)
        {
            std::cout << " -- Solving the system ... " << std::flush;
        }
        linearSolver.solveSystem (rhsBC, solution, systemMatrix);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        // Error computation
        if (verbose)
        {
            std::cout << " -- Computing the error ... " << std::flush;
        }
		//Calcolo errore e norme, abbastanza semplice da leggere
        vector_Type solutionErr (solution);
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( exactSolution ), solutionErr, 0.0 );
        solutionErr -= solution;
        solutionErr.abs();
        Real l2error (uFESpace->l2Error (exactSolution, vector_Type (solution, Repeated), 0.0) );
        if (verbose)
        {
            std::cout << " -- done ! " << std::endl;
        }
        if (verbose)
        {
            std::cout << " ---> Norm L2  : " << l2error << std::endl;
        }
        Real linferror ( solutionErr.normInf() );
        if (verbose)
        {
            std::cout << " ---> Norm Inf : " << linferror << std::endl;
        }

        // Exporter definition and use
        if (verbose)
        {
            std::cout << " -- Defining the exporter ... " << std::flush;
        }

#ifdef HAVE_HDF5
        ExporterHDF5<mesh_Type> exporter ( dataFile, "solution" );
#else
        ExporterVTK<mesh_Type> exporter ( dataFile, "solution" );
#endif
        exporter.setMeshProcId ( meshPtr, Comm->MyPID() ) ;
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Defining the exported quantities ... " << std::flush;
        }
        boost::shared_ptr<vector_Type> solutionPtr (new vector_Type (solution, Repeated) );
        boost::shared_ptr<vector_Type> solutionErrPtr (new vector_Type (solutionErr, Repeated) );
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Updating the exporter ... " << std::flush;
        }
        exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "solution", uFESpace, solutionPtr, UInt (0) );
        exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "error", uFESpace, solutionErrPtr, UInt (0) );
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if (verbose)
        {
            std::cout << " -- Exporting ... " << std::flush;
        }
        exporter.postProcess (0);
        if (verbose)
        {
            std::cout << " done ! " << std::endl;
        }

        if ( std::fabs ( l2error - 0.0006169843149652788l ) > 1e-10 )
        {
            std::cout << " <!> Solution has changed !!! <!> " << std::endl;
            return EXIT_FAILURE;
        }
        if ( std::fabs ( linferror - 0.0001092814405985187l ) > 1e-10 )
        {
            std::cout << " <!> Solution has changed !!! <!> " << std::endl;
            return EXIT_FAILURE;
        }

        if (verbose)
        {
            std::cout << "End Result: TEST PASSED" << std::endl;
        }

    } // needed to properly destroy all objects inside before mpi finalize

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}
