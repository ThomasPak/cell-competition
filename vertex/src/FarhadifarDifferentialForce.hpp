#ifndef FARHADIFARDIFFERENTIALFORCE_HPP_
#define FARHADIFARDIFFERENTIALFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellLabel.hpp"

#include <iostream>

/**
 * A force class for use in Vertex-based simulations. This force is based on the
 * Energy function proposed by Farhadifar et al in  Curr. Biol., 2007, 17, 2095-2104.
 */


template<unsigned DIM>
class FarhadifarDifferentialForce : public AbstractForce<DIM>
{
friend class TestForces;

private:

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mCellAreaElasticityParameter;
        archive & mLabelledCellAreaElasticityParameter;
        archive & mCellPerimeterContractilityParameter;
        archive & mLabelledCellPerimeterContractilityParameter;
        archive & mCellCellLineTensionParameter;
        archive & mLabelledCellLabelledCellLineTensionParameter;
        archive & mLabelledCellCellLineTensionParameter;
        archive & mCellBoundaryLineTensionParameter;
        archive & mLabelledCellBoundaryLineTensionParameter;
    }

protected:

    /**
     * The strength of the area term in the model. Corresponds to K_alpha in
     * Farhadifar's paper.
     */
    double mCellAreaElasticityParameter;

    /**
     * The strength of the area term in the model for a labelled cell.
     * */
    double mLabelledCellAreaElasticityParameter;

    /**
     * The strength of the perimeter term in the model. Corresponds to
     * Gamma_alpha in Farhadifar's paper.
     */
    double mCellPerimeterContractilityParameter;

    /**
     * The strength of the perimeter term in the model for a labelled cell.
     */
    double mLabelledCellPerimeterContractilityParameter;

    /**
     * The strength of the cell-cell line tension between two cells.
     * Lambda_{i,j} in Farhadifar's paper.
     */
    double mCellCellLineTensionParameter;

    /**
     * The strength of the cell-cell line tension between two labelled cells.
     */
    double mLabelledCellLabelledCellLineTensionParameter;

    /**
     * The strength of the cell-cell line tension between a labelled and
     * non-labelled cell.
     */
    double mLabelledCellCellLineTensionParameter;

    /**
     * The strength of the line tension at the boundary. This term does
     * correspond to Lambda_{i,j} in Farhadifar's paper.
     */
    double mCellBoundaryLineTensionParameter;

    /**
     * The strength of the line tension at the boundary for labelled cells.
     */
    double mLabelledCellBoundaryLineTensionParameter;

public:

    /**
     * Constructor.
     */
    FarhadifarDifferentialForce();

    /**
     * Destructor.
     */
    virtual ~FarhadifarDifferentialForce() override;

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the energy function
     * Farhadifar's model.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation) override;

    /**
     * Get the line tension parameter for the edge between two given nodes.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     * @param rVertexCellPopulation reference to the cell population
     *
     * @return the line tension parameter for this edge.
     */
    double GetLineTensionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation);

    /**
     * @return mCellAreaElasticityParameter
     */
    double GetCellAreaElasticityParameter();

    /**
     * @return mLabelledCellAreaElasticityParameter
     */
    double GetLabelledCellAreaElasticityParameter();

    /**
     * @return mCellPerimeterContractilityParameter
     */
    double GetCellPerimeterContractilityParameter();

    /**
     * @return mLabelledCellPerimeterContractilityParameter
     */
    double GetLabelledCellPerimeterContractilityParameter();

    /**
     * @return mCellCellLineTensionParameter
     */
    double GetCellCellLineTensionParameter();

    /**
     * @return mLabelledCellLabelledCellLineTensionParameter
     */
    double GetLabelledCellLabelledCellLineTensionParameter();

    /**
     * @return mLabelledCellCellLineTensionParameter
     */
    double GetLabelledCellCellLineTensionParameter();

    /**
     * @return mCellBoundaryLineTensionParameter
     */
    double GetCellBoundaryLineTensionParameter();

    /**
     * @return mLabelledCellBoundaryLineTensionParameter
     */
    double GetLabelledCellBoundaryLineTensionParameter();

    /**
     * Set mCellAreaElasticityParameter.
     *
     * @param cellAreaElasticityParameter the new value of
     * mCellAreaElasticityParameter
     */
    void SetCellAreaElasticityParameter(double cellAreaElasticityParameter);

    /**
     * Set mLabelledCellAreaElasticityParameter.
     *
     * @param labelledCellAreaElasticityParameter the new value of
     * mLabelledCellAreaElasticityParameter
     */
    void SetLabelledCellAreaElasticityParameter(double labelledCellAreaElasticityParameter);

    /**
     * Set mCellPerimeterContractilityParameter.
     *
     * @param cellPerimeterContractilityParameter the new value of
     * cellPerimeterContractilityParameter
     */
    void SetCellPerimeterContractilityParameter(double cellPerimeterContractilityParameter);

    /**
     * Set mLabelledCellPerimeterContractilityParameter.
     *
     * @param labelledCellPerimeterContractilityParameter the new value of
     * labelledCellPerimeterContractilityParameter
     */
    void SetLabelledCellPerimeterContractilityParameter(
            double labelledCellPerimeterContractilityParameter);

    /**
     * Set mCellCellLineTensionParameter.
     *
     * @param cellLineTensionParameter the new value of mCellCellLineTensionParameter
     */
    void SetCellCellLineTensionParameter(double cellCellLineTensionParameter);

    /**
     * Set mLabelledCellLabelledCellLineTensionParameter.
     *
     * @param labelledCellLabelledCellLineTensionParameter the new value of
     * mLabelledCellLabelledCellLineTensionParameter
     */
    void SetLabelledCellLabelledCellLineTensionParameter(
            double labelledCellLabelledCellLineTensionParameter);

    /**
     * Set mLabelledCellCellLineTensionParameter.
     *
     * @param labelledCellCellLineTensionParameter the new value of
     * mLabelledCellCellLineTensionParameter
     */
    void SetLabelledCellCellLineTensionParameter(
            double labelledCellCellLineTensionParameter);

    /**
     * Set mCellBoundaryLineTensionParameter.
     *
     * @param cellBoundaryLineTensionParameter the new value of mCellBoundaryLineTensionParameter
     */
    void SetCellBoundaryLineTensionParameter(double cellBoundaryLineTensionParameter);

    /**
     * Set mLabelledCellBoundaryLineTensionParameter.
     *
     * @param labelledCellBoundaryLineTensionParameter the new value of
     * mLabelledCellBoundaryLineTensionParameter
     */
    void SetLabelledCellBoundaryLineTensionParameter(
            double labelledCellBoundaryLineTensionParameter);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FarhadifarDifferentialForce)

#endif /*FARHADIFARDIFFERENTIALFORCE_HPP_*/
