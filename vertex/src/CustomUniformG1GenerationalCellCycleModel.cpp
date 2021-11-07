#include "CustomUniformG1GenerationalCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

CustomUniformG1GenerationalCellCycleModel::CustomUniformG1GenerationalCellCycleModel()
    : AbstractSimpleGenerationalCellCycleModel(),
      mRange(2.0)
{
}

CustomUniformG1GenerationalCellCycleModel::CustomUniformG1GenerationalCellCycleModel(const CustomUniformG1GenerationalCellCycleModel& rModel)
   : AbstractSimpleGenerationalCellCycleModel(rModel),
      mRange(rModel.mRange)
{
    /*
     * The member variables mGeneration and mMaxTransitGeneration are
     * initialized in the AbstractSimpleGenerationalCellCycleModel
     * constructor.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* CustomUniformG1GenerationalCellCycleModel::CreateCellCycleModel()
{
    return new CustomUniformG1GenerationalCellCycleModel(*this);
}

double CustomUniformG1GenerationalCellCycleModel::GetRange()
{
    return mRange;
}

void CustomUniformG1GenerationalCellCycleModel::SetRange(double range)
{
    mRange = range;
}

void CustomUniformG1GenerationalCellCycleModel::SetG1Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    assert(mpCell != nullptr);

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        double g1_duration = GetStemCellG1Duration();
        assert(g1_duration - 0.5 * GetRange() >= 0);

        double rand_num = p_gen->ranf();
        mG1Duration = g1_duration - 0.5 * GetRange() + GetRange() * rand_num; // U[13,15] for default parameters (mStemCellG1Duration)
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        double g1_duration = GetTransitCellG1Duration();
        assert(g1_duration - 0.5 * GetRange() >= 0);

        double rand_num = p_gen->ranf();
        mG1Duration = g1_duration - 0.5 * GetRange() + GetRange() * rand_num; // U[1,3] for default parameters (mStemCellG1Duration)
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

void CustomUniformG1GenerationalCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
     *rParamsFile << "\t\t\t<Range>" << mRange << "</Range>\n";

    // No new parameters to output, so just call method on direct parent class
    AbstractSimpleGenerationalCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(CustomUniformG1GenerationalCellCycleModel)
