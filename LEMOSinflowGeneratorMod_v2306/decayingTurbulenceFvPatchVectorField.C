/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) held by original author
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

\*---------------------------------------------------------------------------*/

#include "decayingTurbulenceFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
//#include "ListListOps.H"
#include "PstreamReduceOps.H"
#include "OFstream.H"
#include "IFstream.H"

#include <ctime>



Foam::Random ranGen(Foam::label(time(0)));

const Foam::scalar OVERLAP = 0.3;

const Foam::label NUM_VORT = 10;

inline Foam::scalar lspot(Foam::scalar l)
{
    return 3*l;
}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const fvPatch& p, 
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    LField_      (p.size()),
    refField_    (p.size()),
    RField_      (p.size()),
    curTimeIndex_(-1),
    vortons_     (),
    Lund_        (),
    R_           (),
    ind_         (2),
    direction_   (1),
    mapperPtr_(),
    perturb_(0),
    nVortField_(p.size(),Zero)
{
    //Info << "NULL constructor " << this->dimensionedInternalField().name() << endl;
}


Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    LField_      (ptf.LField_, mapper),
    refField_    (ptf.refField_, mapper),
    RField_      (ptf.RField_, mapper),
    curTimeIndex_(-1),
    vortons_     (ptf.vortons_),
    Lund_        (ptf.Lund_, mapper),
    R_           (ptf.R_, mapper),
    ind_         (ptf.ind_),
    direction_   (ptf.direction_),
    mapperPtr_(),
    perturb_(ptf.perturb_),
    nVortField_  (ptf.nVortField_, mapper)
{

    //Info << "mapper constructor " << this->dimensionedInternalField().name() << endl;
    //Info << "mapper constructor (no field info)"<< endl;
}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const fvPatch& p, 
    const DimensionedField<vector, volMesh>& iF, 
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF),
    LField_      (p.size(),pTraits<scalar>::zero),
    refField_    (p.size(),pTraits<vector>::zero),
    RField_      (p.size(),pTraits<symmTensor>::zero),
    curTimeIndex_(-1),
    Lund_        (p.size(), pTraits<tensor>::zero),
    R_           (p.size(), pTraits<symmTensor>::zero),
    ind_         (2),
    direction_   (readScalar(dict.lookup("direction"))),
    mapperPtr_(),
    perturb_(dict.lookupOrDefault("perturb", 1e-5)),
    nVortField_  (p.size(),pTraits<scalar>::zero)
{
    //Info << "dict constructor " << this->dimensionedInternalField().name() << endl;
    if (dict.found("vortons"))
        vortons_ = SLList<decayingVorton>(dict.lookup("vortons"));

    if (Pstream::master())
    {
        fileName inputFile(db().time().path()/db().time().timeName()/+"vortons_"+this->internalField().name());
        IFstream os2(inputFile);

        if (os2.opened())
        {
            word name;
            os2.read(name);
            if (name=="vortons")
            {
                Info << "found vortons file: "<<name << endl;
                os2 >> vortons_;
            }
        }
    }
    
    //Fills  LField_, refField_, RField_
    readInletData();

    // debug: write the interpolated source field if required.
    if (dict.lookupOrDefault("writeSourceFields",false)==true)
    {
        writeSourceFields();
    }

    R_= RField_;
    if (dict.found("R"))
        R_ = Field<symmTensor>("R", dict, p.size());

    if (dict.found("value"))
        fixedValueFvPatchField<vector>::operator==(Field<vector>("value", dict, p.size()));
    else
        fixedValueFvPatchField<vector>::operator==(refField_);

    forAll(RField_, I)
    {
        if ((RField_[I].xx() < 0) || (RField_[I].yy() < 0) || (RField_[I].zz() < 0))
        FatalErrorIn("decayingTurbulenceFvPatchVectorField::")<< "Some RMS are negative in R"<<abort(FatalError);
    }
    
    bool bError = false;
    bool bRobustLund = dict.lookupOrDefault("robustLundTransformation",false);

    forAll(RField_,i)
    {
        if (RField_[i].component(symmTensor::XX)<=0.0)
        {
            Info << "Found a negative or 0 Re Stress in XX" << endl;
            bError=true;
            Lund_[i]=pTraits<tensor>::I;
        }
    }
    if (!bRobustLund and !bError)
    {
        Info << "Standard Lund transformation algorithm applied." << endl;
        Lund_.replace(tensor::XX, sqrt(RField_.component(symmTensor::XX)));
        Lund_.replace(tensor::YX, RField_.component(symmTensor::XY)/Lund_.component(tensor::XX));
        Lund_.replace(tensor::ZX, RField_.component(symmTensor::XZ)/Lund_.component(tensor::XX));
        Lund_.replace(tensor::YY, sqrt(RField_.component(symmTensor::YY)-sqr(Lund_.component(tensor::YX))));
        Lund_.replace(tensor::ZY, (RField_.component(symmTensor::YZ) - Lund_.component(tensor::YX)*Lund_.component(tensor::ZX) )/Lund_.component(tensor::YY));
        Lund_.replace(tensor::ZZ, sqrt(RField_.component(symmTensor::ZZ) - sqr(Lund_.component(tensor::ZX))-sqr(Lund_.component(tensor::ZY))));
    }
    else if (bRobustLund and !bError)
    {
        Info << "Robust Lund transformation algorithm applied." << endl;
        forAll(Lund_,i)
        {
            bool errorFound = false;
            scalar a11 = sqrt(RField_[i][0]);
            scalar a21 = RField_[i][1]/a11;
            scalar a31 = RField_[i][2]/a11;

            if (RField_[i][3]<=sqr(a21))  // test if a22 exists (risk of sqrt a negative number)
            {
                errorFound = true;
            };
            if (errorFound==false) // test if a32 exists (a22 must be different than zero)
            {
                scalar a22 = sqrt(RField_[i][3]-sqr(a21));
                if (a22==0.0)
                {
                    errorFound = true;
                };
            };
            if (errorFound==false) // test if a33 exists (risk of sqrt a negative number)
            {
                scalar a22 = sqrt(RField_[i][3]-sqr(a21));
                scalar a32 = (RField_[i][4]-a21*a31)/a22;
                if (RField_[i][5]<(sqr(a31)+sqr(a32)))
                {
                    errorFound = true;
                };
            };

            if (errorFound==true)
            {
                Lund_[i][0] = sqrt(RField_[i][0]);  //xx
                Lund_[i][4] = sqrt(RField_[i][3]);  //yy
                Lund_[i][8] = sqrt(RField_[i][5]);  //zz
            }
            else
            {
                Lund_[i][0] = sqrt(RField_[i][0]);  //xx
                Lund_[i][3] = RField_[i][1]/Lund_[i][0];  //yx
                Lund_[i][4] = sqrt(RField_[i][3]-sqr(Lund_[i][3]));  //yy
                Lund_[i][6] = RField_[i][2]/Lund_[i][0];  //zx
                Lund_[i][7] = (RField_[i][4]-Lund_[i][3]*Lund_[i][6])/Lund_[i][4];  //zy
                Lund_[i][8] = sqrt(RField_[i][5]-sqr(Lund_[i][6])-sqr(Lund_[i][7]));  //zz
            };

        };
    }else
    {
        Lund_=pTraits<tensor>::I;
    }

    // debug: write the interpolated source field if required.
    if (dict.lookupOrDefault("writeSourceFields",false)==true)
    {
        fileName rootPath(this->db().time().constant()/"boundaryData"/this->patch().name());
        OFstream(rootPath/"LundPatch")() << Lund_;
    }

    if (dict.found("ind"))
        ind_ = readScalar(dict.lookup("ind"));

    const scalarField& sf = patch().magSf();
    
    int numsf = 0;
    
    forAll(sf, I)
    {
        if (sf[I] > LField_[I]*LField_[I])
	    numsf++;
    }
    
    if (numsf > 0)
        Pout<<"Warning: there are "<<numsf<<" inlet cells (out of "<<patch().size()<<") where the integral length is smaller than the cell size."<<endl;
}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    LField_      (ptf.LField_),
    refField_    (ptf.refField_),
    RField_      (ptf.RField_),
    curTimeIndex_(-1),
    vortons_     (ptf.vortons_),
    Lund_        (ptf.Lund_),
    R_           (ptf.R_),
    ind_         (ptf.ind_),
    direction_   (ptf.direction_),
    mapperPtr_(),
    perturb_(ptf.perturb_),
    nVortField_  (ptf.nVortField_)
{
    //Info << "ptf constructor " << this->dimensionedInternalField().name() << endl;
}

Foam::decayingTurbulenceFvPatchVectorField::decayingTurbulenceFvPatchVectorField
(
    const decayingTurbulenceFvPatchVectorField& ptf, 
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    LField_      (ptf.LField_),
    refField_    (ptf.refField_),
    RField_      (ptf.RField_),
    curTimeIndex_(-1),
    vortons_     (ptf.vortons_),
    Lund_        (ptf.Lund_),
    R_           (ptf.R_),
    ind_         (ptf.ind_),
    direction_   (ptf.direction_),
    mapperPtr_(),
    perturb_(ptf.perturb_),
    nVortField_  (ptf.nVortField_)
{

    //Info << "ptf,If constructor " << this->dimensionedInternalField().name() << endl;
}

void Foam::decayingTurbulenceFvPatchVectorField::autoMap(const fvPatchFieldMapper& m)
{
    //Info << "autoMap " << this->dimensionedInternalField().name() << endl;

    Field<vector>::autoMap(m);
    LField_.autoMap(m);
    refField_.autoMap(m);
    RField_.autoMap(m);
    R_.autoMap(m);
    nVortField_.autoMap(m);

    // Clear interpolator
    mapperPtr_.clear();

}

void Foam::decayingTurbulenceFvPatchVectorField::rmap(const fvPatchField<vector>& ptf, const labelList& addr)
{
    //Info << "rmap " << this->dimensionedInternalField().name() << endl;
    //Info << "rmap " << endl;
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const decayingTurbulenceFvPatchVectorField& tiptf = refCast<const decayingTurbulenceFvPatchVectorField>(ptf);

    LField_.rmap(tiptf.LField_, addr);
    refField_.rmap(tiptf.refField_, addr);
    RField_.rmap(tiptf.RField_, addr);
    R_.rmap(tiptf.R_, addr);
    nVortField_.rmap(tiptf.nVortField_, addr);


    // Clear interpolator
    mapperPtr_.clear();

}


void Foam::decayingTurbulenceFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
        return;


    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        doUpdate();
        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}

void Foam::decayingTurbulenceFvPatchVectorField::readInletData()
{
    // Initialise
    if (mapperPtr_.empty())
    {
        pointIOField samplePoints
        (
            IOobject
            (
                "points",
                this->db().time().constant(),
                "boundaryData"/this->patch().name(),
                this->db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                false
            )
        );

        // Allocate the interpolator
        mapperPtr_.reset
        (
            new pointToPointPlanarInterpolation
            (
                samplePoints,
                 this->patch().patch().faceCentres(),
                perturb_
            )
        );

        // Reread values and interpolate
        vectorIOField ref_
        (
            IOobject
            (
                "ref",
                this->db().time().constant(),
                "boundaryData"
               /this->patch().name(),
                this->db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                false
            )
        );

        // Reread values and interpolate
        IOField<scalar> L_
        (
            IOobject
            (
                "L",
                this->db().time().constant(),
                "boundaryData"
               /this->patch().name(),
                this->db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                false
            )
        );

        // Reread values and interpolate
        IOField<symmTensor> Rs_
        (
            IOobject
            (
                "R",
                this->db().time().constant(),
                "boundaryData"
               /this->patch().name(),
                this->db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                false
            )
        );

        LField_ = mapperPtr_().interpolate(L_);
        refField_ = mapperPtr_().interpolate(ref_);
        RField_ = mapperPtr_().interpolate(Rs_);

    }
}

void Foam::decayingTurbulenceFvPatchVectorField::writeSourceFields()
{
    // example to save data in files with header (not used)
//    IOField<symmTensor> IORField
//    (
//        IOobject
//        (
//            "RPatch_withHeader",
//            this->db().time().constant(),
//            "boundaryData"
//           /this->patch().name(),
//            this->db(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE,
//            false
//        ),
//        RField_
//    );
//    IORField.write();

    // save data in files without header
    fileName rootPath(this->db().time().constant()/"boundaryData"/this->patch().name());

    OFstream(rootPath/"pointsPatch")() << this->patch().Cf();
    OFstream(rootPath/"RPatch")() << RField_;
    OFstream(rootPath/"LPatch")() << LField_;
    OFstream(rootPath/"refPatch")() << refField_;
}

void Foam::decayingTurbulenceFvPatchVectorField::doUpdate()
{
    Field<vector>& patchField = *this;
    
    Field<vector> tloc = patch().Cf();
    
    List<List<scalar> > l(Pstream::nProcs());
    List<List<vector> > cf(Pstream::nProcs());
    List<List<vector> > rf(Pstream::nProcs());
    
    l[Pstream::myProcNo()].setSize(LField_.size());
    cf[Pstream::myProcNo()].setSize(tloc.size());
    rf[Pstream::myProcNo()].setSize(refField_.size());
    
    label i = 0;
    forAll(tloc, I)
    {
        l[Pstream::myProcNo()][i]  = LField_[I];
        cf[Pstream::myProcNo()][i] = tloc[I];
        rf[Pstream::myProcNo()][i] = refField_[I];

        ++i;
    }
    
    Pstream::gatherList(l);
    Pstream::gatherList(cf);
    Pstream::gatherList(rf);
    
    List<scalar> L  = ListListOps::combine<List<scalar> >(l, accessOp<List<scalar> >());
    List<vector> CF = ListListOps::combine<List<vector> >(cf, accessOp<List<vector> >());
    List<vector> RF = ListListOps::combine<List<vector> >(rf, accessOp<List<vector> >());
    
    List<List<decayingVorton> > v(Pstream::nProcs());
    
    if (Pstream::master())
    {
        Pout<<"Starting inflgen, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;
        forAll(CF, I)
        {
            if (L[I] > 0.0)
            {
                scalar x    = CF[I].x();
                scalar ls   = lspot(L[I]);
                scalar ymin = CF[I].y() - L[I];
                scalar zmin = CF[I].z() - L[I];

                for (label k = 0; k < NUM_VORT; k++)
                {
                    vector v((direction_ > 0) ? x-ls : x+ls, ymin+ranGen.sample01<scalar>()*2*L[I], zmin+ranGen.sample01<scalar>()*2*L[I]);
                    bool allowed = true;
                    for (SLList<decayingVorton>::const_iterator it = vortons_.begin(); it != vortons_.end(); ++it)
                    {
                        if (mag(v - it().location()) < OVERLAP*ls)
                        {
                            allowed = false;
                            break;
                        }
                    }
                    if (!allowed)
                        continue;
                    else
                    {
                        vortons_.insert(decayingVorton(L[I], v, RF[I], (direction_ > 0) ? x+ls : x-ls));
                        if(debug)
                        {
                            nVortField_[I]++;
                        }
                    }
                }
            }
        }

        Pout<<"The number of vortons is "<<vortons_.size()<<endl;

        v[Pstream::myProcNo()].setSize(vortons_.size());
        i = 0;
        for (SLList<decayingVorton>::iterator it = vortons_.begin(); it != vortons_.end(); ++it)
            v[Pstream::myProcNo()][i++] = it();
    }
    
    Pstream::scatterList(v);
    
    List<decayingVorton> V = ListListOps::combine<List<decayingVorton> >(v, accessOp<List<decayingVorton> >());
    
    if (size() > 0)
    {
        Field<vector> turbulent(refField_.size(), pTraits<vector>::zero);
        
        forAll(patchField, I)
        {
            vector pt = tloc[I];
            
            forAll(V, vI)
            {
                const decayingVorton& vorton = V[vI];

                if (mag(vorton.location() - pt) < lspot(vorton.length())) 
                    turbulent[I] += vorton.velocityAt(pt);
            }
        }

        R_ = ((ind_ - 1.0)/ind_)*R_ + (1.0/ind_)*sqr(turbulent);

        Field<tensor> C_(R_.size(), pTraits<tensor>::zero);
        forAll(C_, I)
        {
            C_[I].xx() = 1.0;
            C_[I].yy() = 1.0;
            C_[I].zz() = 1.0;
        }

        forAll(R_, I)
        {
            scalar D1 = R_[I].xx();
            if (D1 > 0)
                C_[I].xx() = 1.0/sqrt(D1);
    
            scalar D2 = R_[I].xx()*R_[I].yy() - R_[I].xy()*R_[I].xy();
            if (D1 > 0 && D2 > 0)
            {
                C_[I].yx() = -R_[I].xy()/sqrt(D1*D2);
                C_[I].yy() = sqrt(D1/D2);
            }

            scalar D3 = det(R_[I]);
            if (D2 > 0 && D3 > 0)
            {
                C_[I].zx() = (R_[I].xy()*R_[I].yz()-R_[I].yy()*R_[I].xz())/sqrt(D2*D3);
                C_[I].zy() = -(R_[I].xx()*R_[I].yz()-R_[I].xz()*R_[I].xy())/sqrt(D2*D3);
                C_[I].zz() = sqrt(D2/D3);
            }
        }

        ind_++;
        
        turbulent = C_&turbulent;
        turbulent = Lund_&turbulent;

        fixedValueFvPatchField<vector>::operator==(refField_+ turbulent);
    }
    
    if (Pstream::master())
    {
        for (SLList<decayingVorton>::iterator it = vortons_.begin(); it != vortons_.end(); ++it)
	    it().move(this->db().time().deltaT().value());
	
        bool modified;
        do
        {
            modified=false;
            for (SLList<decayingVorton>::iterator it = vortons_.begin(); it!=vortons_.end(); ++it)
            {   
	        if (direction_ > 0)
		{
                    if (it().location().x() > it().xmax())
                    {
                        vortons_.remove(it);
                        modified=true;
                        break;
                    }
		}
		else
		{
                    if (it().location().x() < it().xmax())
                    {
                        vortons_.remove(it);
                        modified=true;
                        break;
                    }
		}
            } 
        } while (modified);
	
        Pout<<"Finishing inflgen, time = "<<this->db().time().elapsedClockTime()<<" s"<<endl;
    }
}

void Foam::decayingTurbulenceFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    //LField_.writeEntry("LField", os);
    //refField_.writeEntry("refField", os);
    //RField_.writeEntry("RField", os);
    os.writeKeyword("direction")<<direction_<<token::END_STATEMENT<<nl;
    this->writeEntry("value", os);
    R_.writeEntry("R", os);
    os.writeKeyword("ind")<<ind_<<token::END_STATEMENT<<nl;

    if (Pstream::master())
    {
        fileName outputFile(db().time().path()/db().time().timeName()/+"vortons_"+this->internalField().name());
        OFstream os2(outputFile);
        os2.writeKeyword("vortons")<<vortons_<<token::END_STATEMENT<<nl;

        if(debug)
        {
            fileName rootPath(this->db().time().timePath());
            OFstream(rootPath/"nVortField")() << nVortField_;

        }
    }

    //if (Pstream::master())
    //   os.writeKeyword("vortons")<<vortons_<<token::END_STATEMENT<<nl;
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        decayingTurbulenceFvPatchVectorField
    );
}

