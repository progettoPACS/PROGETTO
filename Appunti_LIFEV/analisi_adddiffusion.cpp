/*
	ANALISI DEI METODI DELLA CLASSE ADRASSEMBLER 
*/

//	DIFFUSIONE
    /*
      This method adds the diffusion matrix times the coefficient
      (whose default value is 1.0) to the matrix passed in argument.
      @Remark The matrix is NOT finalized, you have to call globalAssemble
      outside this method when the matrix is finished.
     */

	/*
		Questa è la funzione che viene chiamata nel main
	*/
    inline void addDiffusion (matrix_ptrType matrix, const Real& coefficient = 1.0)
    {
        addDiffusion (matrix, coefficient, 0, 0);
    }

    /*
	Questa è la chiamata vera che viene eseguita con termini aggiuntivi.
      This method adds the diffusion in the given matrix. Additional arguements are provided
      so that one can chose where to add the mass in the matrix, i.e. somewhere else
      than in the upper left block.
     */
    void addDiffusion (matrix_ptrType matrix, const Real& coefficient, const UInt& offsetLeft, const UInt& offsetUp);

template<typename mesh_type, typename matrix_type, typename vector_type>
void
ADRAssembler< mesh_type, matrix_type, vector_type>::
addDiffusion (matrix_ptrType matrix, const Real& coefficient, const UInt& offsetLeft, const UInt& offsetUp)
{
    // Check that the fespace is set
    ASSERT (M_fespace != 0, "No FE space for assembling the diffusion!");

    M_diffusionAssemblyChrono.start();

    // Some constants
    const UInt nbElements (M_fespace->mesh()->numElements() );	//Elementi della mesh
    const UInt fieldDim (M_fespace->fieldDim() );		//Dimensione del dominio (1,2,3)
    const UInt nbTotalDof (M_fespace->dof().numTotalDof() );	//Numero totale di gradi di libertà

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        //CONTO DELLE QUANTITÀ CHE MI SERVONO PER FARE GLI INTEGRALI APPOGGIANDOMI SULLA CLASSE CurrentFE

        M_diffCFE->update ( M_fespace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );
	/* 
		M_diffCFE: variabile del private, (CurrentFE for the diffusion)
		currentFE_ptrType   boost::scoped_ptr<currentFE_type>  CurrentFE
		CurrentFE è una classe
	
		In pratica sto chiamando il metodo update della classe CurrentFE sull'elementino di mesh con le flag che
		aggiornano i valori delle derivate delle basi e i pesi di quadratura (per la diffusione serve solo quello)
	*/
        // Clean the local matrix
        M_localDiff->zero();
	//Non è altro che un puntatore alla matrice locale che utilizzo per costruire la matrice di stiffness
	//Immagino che il metodo zero crei la matrice piena di zeri.

        // local stiffness
        AssemblyElemental::stiffness (*M_localDiff, *M_diffCFE, coefficient, fieldDim);

//! Elementary stiffness for constant coefficient
/*!
  This function assembles the local stiffness matrix when the coefficient is constant.

  @param localStiff The local matrix to be filled (not cleaned by this function)
  @param stiffCFE The currentFE structure already updated for the assembly. It requires
  dphi and wDetJacobian to be accessible.
  @param coefficient The coefficient
  @param fieldDim The dimension of the FE space (scalar/vectorial)
 */



        /*
		 Itera sulle dimensioni del modello quindi in realtà si crea una matricetta per ogni dimensione,
		e con il metodo assembleMatrix le va a posizionare all'interno del matricione passato alla funzione
		addDiffusion.
	*/
        for (UInt iFieldDim (0); iFieldDim < fieldDim; ++iFieldDim)
        {
            assembleMatrix ( *matrix,
                             *M_localDiff,
                             *M_diffCFE,
                             *M_diffCFE,
                             M_fespace->dof(),
                             M_fespace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * nbTotalDof + offsetLeft, iFieldDim * nbTotalDof + offsetUp );
        }
 }

    M_diffusionAssemblyChrono.stop();
}


//---------------------------- ANALISI DELLA FUNZIONE STIFFNESS --------------------------------------------


void stiffness (MatrixElemental& localStiff,
                const CurrentFE& stiffCFE,
                const Real& coefficient,
                const UInt& fieldDim)
{
    const UInt nbFEDof (stiffCFE.nbFEDof() );
    const UInt nbQuadPt (stiffCFE.nbQuadPt() );
    Real localValue (0);

    MatrixElemental::matrix_type mat_tmp (nbFEDof, nbFEDof);

    // Loop over the basis functions (diagonal part)
    for (UInt iDof (0); iDof < nbFEDof; ++iDof)
    {
        localValue = 0.0;

        // Loop on the quadrature noodes
        for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
        {
            for (UInt iDim (0); iDim < stiffCFE.nbLocalCoor(); ++iDim)
            {
                localValue += stiffCFE.dphi (iDof, iDim, iQuadPt)
                              * stiffCFE.dphi (iDof, iDim, iQuadPt)
                              * stiffCFE.wDetJacobian (iQuadPt);
/*
	Ciclo con due/tre for annidati per scorrere tutti i nodi del volume della cella
	Questo pezzo è altamente parallelizzabile, viene eseguito in parallelo o no?

	+ stiffCFE_Himod.dphi()
	 *stiffCFE_Himod.dphi()
	 *stiffCFE.phi()
	 *stiffCFE.phi()
	 *stiffCFE_Himod.wDetJacobian()
	 *stiffCFE:wDetJacobian()
        + stiffCFE_Himod.phi()
         *stiffCFE_Himod.phi()
         *stiffCFE.dphi()
         *stiffCFE.dphi()
         *stiffCFE_Himod.wDetJacobian()
         *stiffCFE:wDetJacobian()
*/
            }
        }
        localValue *= coefficient;
/*
	Per quanto riguarda un coefficiente di diffusione dipendente dallo spazio bisognerebbe comportarsi similmente al
	caso delle basi Himod e fare una classe stiffCFE_funzioni in modo tale da inserire nell'integrale pure lei.
*/
        // Add on the local matrix
        mat_tmp (iDof, iDof) = localValue;
    }
 // Loop over the basis functions (extradiagonal part)
    for (UInt iDof (0); iDof < nbFEDof ; ++iDof)
    {
        for (UInt jDof (0); jDof < iDof; ++jDof)
        {
            localValue = 0.0;

            // Loop on the quadrature nodes
            for (UInt iQuadPt (0); iQuadPt < nbQuadPt; ++iQuadPt)
            {
                for (UInt iDim (0); iDim < stiffCFE.nbLocalCoor(); ++iDim)
                {
                    localValue += stiffCFE.dphi (iDof, iDim, iQuadPt)
                                  * stiffCFE.dphi (jDof, iDim, iQuadPt)
                                  * stiffCFE.wDetJacobian (iQuadPt);
                }
            }


            localValue *= coefficient;

            // Put in the local matrix
            mat_tmp (iDof, jDof) = localValue;
            mat_tmp (jDof, iDof) = localValue;
        }
    }

    // Copying the diffusion in all the diagonal blocks
    for ( UInt iDim (0); iDim < fieldDim; ++iDim)
    {
        MatrixElemental::matrix_view mat = localStiff.block (iDim, iDim);
        mat += mat_tmp;
    }
}

