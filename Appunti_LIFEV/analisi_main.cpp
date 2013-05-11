/*
	ANALISI DEL MAIN RIFERITO AL CASO ADR 1D
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
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}

//--------------------- DEFINIZIONE DI FUNZIONI UTILI AL PROBLEMA PER CONDIZIONI DI BORDO E FORZANTI ------------------------

Real exactSolution ( const Real& /* t */, const Real& x, const Real& /*y*/, const Real& /* z */, const ID& /* i */ )
{
    return std::sin ( M_PI * 0.5 * x );
}


Real fRhs ( const Real& /* t */, const Real& x, const Real& /* y */, const Real& /* z */ , const ID& /* i */ )
{
    return 1.25 * M_PI * M_PI * std::sin ( M_PI * 0.5 * x );
}

//-----------------------------------------------------------------------------------------------------------------------



//-------------------------------------- INIZIO DEL MAIN -------------------------------------------------
int
main ( int argc, char* argv[] )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
#endif

    {
        // needed to properly destroy all objects inside before mpi finalize

#ifdef HAVE_MPI
        boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
        ASSERT ( Comm->NumProc() < 2, "The test does not run in parallel." );
#else
        boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

        typedef RegionMesh<LinearLine> mesh_Type;
        typedef MatrixEpetra<Real> matrix_Type;
        typedef VectorEpetra vector_Type;
        typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
        typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

        const bool verbose (Comm->MyPID() == 0);

//--------------------- LETTURA DATI TRAMITE GETPOT -------------------------------------------

        // Read first the data needed
        GetPot dataFile ( "data" );
  

        // Build the mesh
        MeshData meshData (dataFile, "mesh");
        boost::shared_ptr< mesh_Type > meshPtr ( new mesh_Type ( Comm ) );

        // Set up the structured mesh
        regularMesh1D ( *meshPtr, 0,
                        dataFile ( "mesh/n", 20 ),
                        dataFile ( "mesh/verbose", false ),
                        dataFile ( "mesh/length", 1. ),
                        dataFile ( "mesh/origin", 0. ));

//------------------------ COSTRUZIONE DELLO SPAZIO ELEMENTI FINITI ----------------------------------------

        feSpacePtr_Type uFESpace ( new feSpace_Type ( meshPtr, feSegP1, quadRuleSeg1pt, quadRuleNode1pt, 1, Comm ) );
/*
    FESpace ( meshPtr_Type            mesh,	shared_ptr alla mesh
              const ReferenceFE&      refFE,	Reference FE for the velocity
              const QuadratureRule&   Qr,	Quadrature rule for volumic elementary computations
              const QuadratureRule&   bdQr,	Quadrature rule for surface elementary computations
              const Int               fDim,	Dimensione della variabile di campo (scalare/vettoriale)
              const commPtr_Type&     commptr	createMap(commptr) Penso serva all'inizializzazione delle mappa
            );
*/  


//---------------------- COSTRUZIONE DELL'OGGETTO ASSEMBLER E DELLA MATRICE DI SISTEMA ------------------------------------------
       
	 ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
/*
	È possibile costruire l'oggetto inizializzando immediatamente altrimenti così utilizza i costruttori di default
*/
        adrAssembler.setFespace (uFESpace);
/*
	Ho provato a guardare ma c'è una sequela di typedef indecifrabile o meglio lunghissima, quindi semplicemente
	definito l'oggetto adrAssembler è posibile inisializzarlo con il metodo setFespace anche più volte immagino.	
*/
        boost::shared_ptr<matrix_Type> systemMatrix (new matrix_Type ( uFESpace->map() ) );
/*
	DEFINIZIONE DELLA MATRICE DI SISTEMA
	È uno shared_ptr del tipo MatrixEpetra<Real>.

	COSTRUZIONE DI UNA MATRICE QUADRATA O RETTANGOLARE
	@param map Row map. The column map will be defined in MatrixEpetra<DataType>::GlobalAssemble(...,...)
	@param numEntries The average number of entries for each row.
    
	MatrixEpetra ( const MapEpetra& map, Int numEntries = 50 );
	
	Il metodo map() infatti restituisce il puntatore alla variabile private nell'oggetto uFESpace della mappa creata al
	momento della costruzione dell'oggetto con createMap(commptr).
*/


//	Una volta inizializzato l'oggetto ADRAssembler occorre aggiungerci i pezzi che servono (diffusione,reazione e avvezione)
        adrAssembler.addDiffusion (systemMatrix, 1.);
        adrAssembler.addMass (systemMatrix, M_PI * M_PI );

//	In realtà questa funzione non fa altro che richiamare il GlobalAssemble() delle librerie trilinos immagino (potrei
//	sbagliarmi) e rimette insieme i pezzi in modo tale che la matrice di sistema sia pronta all'utilizzo
        systemMatrix->globalAssemble();



//---------------- COSTRUZIONE DEL VETTORE TERMINI NOTI RHS ----------------------------------------------------

        vector_Type rhs (uFESpace->map(), Repeated);
/*
	This constructor uses maps to build the vector
	@param map Map to be used to split the vector between the processors
	@param mapType Specify wether the map is Unique or Repeated (COSA SI INTENDE CON QUESTI DUE TERMINI???)

	VectorEpetra ( const MapEpetra& map, const MapEpetraType& mapType = Unique );
*/

        vector_Type fInterpolated (uFESpace->map(), Repeated);
        
	uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( fRhs ), fInterpolated, 0.0 );
/*
	Sembra essere una funzione utilità dell'oggetto uFESpace che permette data una funzione del tipo fRhs di interpolarla
	sulla mesh e genereare il vettore noto corrispondente. Nella definizione la terza variabile la chiama time e la usa 
	chiamando la funzione 
	nodalValues[iterDof] =  fct (time, interpCFE.quadNode (iterDof, 0), interpCFE.quadNode (iterDof, 1)
                                             , interpCFE.quadNode (iterDof, 2), iDim);
	dove fct è sostanzialmente fRhs, sicuramente deve trattarsi dell'evoluzione temporale, in questo esempio assente,
	quindi fissata all'istante 0.0
*/

        adrAssembler.addMassRhs (rhs, fInterpolated);
        rhs.globalAssemble();
/*
	Immaggino come nel caso precedente aggiunge i pezzi necessari e gli assembla insieme sul vettor RHS
*/



//------------------- CONDIZIONI DI BORDO DEFINIZIONE E APPLICAZIONE ----------------------------------
//Si trovano nella cartella fem

        BCHandler bchandler;
/*
	È l'oggetto che contiene tutte le condizioni di bordo.
	La variabile è un vettore di BCBase e un set di bc_Flag_Type (sono markerID_Type). BCBase è invece un oggetto,
	amico di BCHandler, che contiene tutti i dati della condizione di bordo e un sacco di funzioni utili quali getter e setter
	per maneggiarla e parlarci.
*/
        BCFunctionBase BCu ( static_cast<feSpace_Type::function_Type> ( exactSolution ) );
/*
	È l'oggetto che tiene la funzione da applicare ai bordi, contiene numerose funzioni utili e operatori, è possibile
	anche risettarla. L'oggetto nelle private è boss::function<Real (Real,Real,Real,Real,ID)>.
*/

        bchandler.addBC ("Dirichlet", Structured1DLabel::LEFT, Essential, Full, BCu, 1);
        bchandler.addBC ("Dirichlet", Structured1DLabel::RIGHT, Essential, Full, BCu, 1);
/*
	Boundary conditions are added using the method 

	void addBC ( const std::string& name,
           	 const markerID_Type& f,
		 const bcType_Type& type,
                 const bcMode_Type& mode,
                 BCFunctionBase& bcFunction,
                 const UInt& numberOfComponents); 

enum bcType_Type
{
    Natural,            < Neumann boundary conditions 
    Robin,              < Robin boundary conditions 
    Flux,               < Flux boundary conditions
    Resistance,         < Resistance boundary conditions
    Essential,          < Dirichlet boundary conditions
    EssentialEdges,     < Dirichlet boundary conditions on edges
    EssentialVertices   < Dirichlet boundary conditions on vertices 
};
enum bcMode_Type
{
    Scalar,     < To be used for scalar problems 
    Full,       < To be used for vector problems, when the boundary condition involves all components
    Component,  < To be used for vector problems, when the boundary condition DOESN'T involve all components
    Normal,     < To be used for vector problems, when the boundary condition involve the normal component
    Tangential, < To be used for vector problems, when the boundary condition involve the tangential component
    Directional < To be used for vector problems, when the boundary condition involve a specific direction
};

	@c type is a @c bcType_Type and describes the type of boundary condition (@c Essential, @c Natural, etc.).
	@c flag is used to determine on which elements the boundary condition should be prescribed.
*/


        bchandler.bcUpdate (*uFESpace->mesh(), uFESpace->feBd(), uFESpace->dof() );
/* 
	Update all the boundary conditions
    
      This method update the BC classes checking the markers
      on the boundary.
      Except for EssentialEdges and EssentialVertices type of boundary conditions,
      the flags associated with the boundary conditions will be searched among the markers of the boundary element (facets).

      Then, all the DOFs belonging to a matching boundary element will be associated with the appropriate boundary conditions.
      It is possible then the same DOF is shared by different boundary conditions.
      In the case of essential boundary conditions, the largest condition will be prescribed on the shared DOF (largest in the ordering given by the operator< in BCBase.hpp)
      In particular, if two Essential boundary conditions share the same DOF, it will be prescribed the condition with the largest flag.
      This behavior is due to the fact that the largest boundary condition is the last to be prescribed.

      Finally M_bcUpdateDone is set to true, and it is possible to prescribed boundary conditions using functions in BCManage.hpp.

      @param mesh The mesh
      @param boundaryFE Current finite element on the boundary
      @param dof Container of the local to global map of DOF id
*/



        vector_Type rhsBC (rhs, Unique);
//	Sarà fatto in modo che aggiunge componenti e basta (piuttosto che sostituire)
        bcManage (*systemMatrix, rhsBC, *uFESpace->mesh(), uFESpace->dof(), bchandler, uFESpace->feBd(), 1.0, 0.0);
        rhs = rhsBC;

/*
	Prescribe boundary conditions.

	Insomma mette insieme le condizione di bordo settate con handler andando a modificare il termine noto e la matrice
	di sistema che viene infatti passata.

  @param matrix   The system matrix
  @param rightHandSide   The system right hand side
  @param mesh  The mesh
  @param dof  Container of the local to global map of DOFs
  @param bcHandler The boundary conditions handler
  @param currentBdFE Current finite element on boundary
  @parma diagonalizeCoef The coefficient used during the system diagonalization
  @param time The time
 */





//---------------------- DEFINIZIONE DEL SOLVER -----------------------------------------------------

        SolverAztecOO linearSolver;
/*
	Ricordarsi che Aztec00 è deprecato quindi non dovremmo utilizzare lui, inoltre se non utilizziamo lui
	bisogna utilizzare un diverso input da file che possiamo vedere nell'esempio su Darcy
*/
        linearSolver.setDataFromGetPot (dataFile, "solver");
        linearSolver.setupPreconditioner (dataFile, "prec");

        linearSolver.setMatrix (*systemMatrix);

        linearSolver.setCommunicator (Comm);

        // Definition of the solution

        vector_Type solution (uFESpace->map(), Unique);

        // Solve the solution

        linearSolver.solveSystem (rhsBC, solution, systemMatrix);


        // Error computation

        vector_Type solutionErr (solution);
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( exactSolution ), solutionErr, 0.0 );
        solutionErr -= solution;
        solutionErr.abs();
        Real l2error (uFESpace->l2Error (exactSolution, vector_Type (solution, Repeated), 0.0) );

        Real linferror ( solutionErr.normInf() );

//------------------------- DEFINIZIONE E UTILIZZO DELL'EXPORT --------------------------------------------

#ifdef HAVE_HDF5
        ExporterHDF5<mesh_Type> exporter ( dataFile, "solution" );
#else
        ExporterVTK<mesh_Type> exporter ( dataFile, "solution" );
#endif
        exporter.setMeshProcId ( meshPtr, Comm->MyPID() ) ;

        boost::shared_ptr<vector_Type> solutionPtr (new vector_Type (solution, Repeated) );
        boost::shared_ptr<vector_Type> solutionErrPtr (new vector_Type (solutionErr, Repeated) );

        exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "solution", uFESpace, solutionPtr, UInt (0) );
        exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "error", uFESpace, solutionErrPtr, UInt (0) );

        exporter.postProcess (0);

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



    } // needed to properly destroy all objects inside before mpi finalize

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}
