/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::decayingTurbulenceFvPatchVectorField

Description
    This boundary condition provides a decayingTurbulence condition,
    calculated as:

        \f[
            Q = Q_{0} + Q_{p} + s*Q_{t}
        \f]

    where
    \vartable
        s      | single scalar value [units]
        Q_{0}  | single vector value [units]
        Q_{p}  | vector field across patch [units]
        Q_{t}  | vector function of time [units]
    \endtable

Usage
    \table
        Property   | Description                  | Req'd? | Default
        scalarData | single scalar value          | yes    |
        data       | single vector value          | yes    |
        fieldData  | vector field across patch    | yes    |
        timeVsData | vector function of time      | yes    |
        wordData   | word, eg name of data object | no     | wordDefault
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type       decayingTurbulence;
        scalarData -1;
        data       (1 0 0);
        fieldData  uniform (3 0 0);
        timeVsData table (
                             (0 (0 0 0))
                             (1 (2 0 0))
                          );
        wordName   anotherName;
        value      uniform (4 0 0); // optional initial value
    }
    \endverbatim

SourceFiles
    decayingTurbulenceFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef BC_H
#define BC_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

#include "Random.H"
#include "SLList.H"
#include "pointToPointPlanarInterpolation.H"

#include "decayingVorton.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class decayingTurbulenceFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class decayingTurbulenceFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{


















    // Private Data

        //- Single valued scalar quantity, e.g. a coefficient
        //scalar scalarData_;

	//- Fraction of perturbation (fraction of bounding box) to add
	scalar perturb_;

        //- Single valued Type quantity, e.g. reference pressure pRefValue_
        //  Other options include vector, tensor
        //vector data_;

        //- Field of Types, typically defined across patch faces
        //  e.g. total pressure p0_.  Other options include vectorField
        //vectorField fieldData_;

	scalarField nVortField_; //debug field, counts number of vortons

	scalarField LField_;

	vectorField refField_;

	tensorField Lund_;

	symmTensorField R_;

	symmTensorField RField_;



        //- Type specified as a function of time for time-varying BCs
        //autoPtr<Function1<vector>> timeVsData_;

	//- 2D interpolation
	autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- Word entry, e.g. pName_ for name of the pressure field on database
        //word wordData_;

        //- Label, e.g. patch index, current time index
        //label labelData_;

        //- Current time index (used for updating)
        label curTimeIndex_;

	label ind_;

	label direction_;

        //- Boolean for true/false, e.g. specify if flow rate is volumetric_
        //bool boolData_;

	SLList<decayingVorton> vortons_;




    // Private Member Functions

        //- Return current time
        scalar t() const;

	void doUpdate();

	void readInletData();

	void writeSourceFields();


public:

    //- Runtime type information
    TypeName("decayingTurbulence");


    // Constructors

        //- Construct from patch and internal field
        decayingTurbulenceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        decayingTurbulenceFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        decayingTurbulenceFvPatchVectorField
        (
            const decayingTurbulenceFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        decayingTurbulenceFvPatchVectorField
        (
            const decayingTurbulenceFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        decayingTurbulenceFvPatchVectorField
        (
            const decayingTurbulenceFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new decayingTurbulenceFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchVectorField&, const labelList&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchVectorField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
