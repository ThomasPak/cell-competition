#ifndef T2SWAPANCESTORISAPOPTOTICWRITER_HPP
#define T2SWAPANCESTORISAPOPTOTICWRITER_HPP

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "CustomVertexBasedCellPopulation.hpp"

/**
 * A writer class to output the time, ancestors and whether the cell was
 * apoptotic of T2 swaps to a file.
 *
 * The output file is called T2AncestorsIsApoptotics.dat by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class T2SwapAncestorIsApoptoticWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    T2SwapAncestorIsApoptoticWriter();

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is defined for use with a CustomVertexBasedCellPopulation only.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) override;

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is defined for use with a CustomVertexBasedCellPopulation only.
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation) override;

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is defined for use with a CustomVertexBasedCellPopulation only.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation) override;

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is defined for use with a CustomVertexBasedCellPopulation only.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation) override;

    /**
     * Visit the population and write the data.
     *
     * This function will attempt to dynamically cast the pointer to
     * CustomVertexBasedCellPopulation*.  If successful, it will use the
     * specialisation function below.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation) override;

    /**
     * Visit the CustomVertexBasedCellPopulation and write the location of any T2 swaps at the present
     * simulation time.
     *
     * Outputs a line of tab-separated values of the form:
     * [num T2 swaps] [T2 swap 0 ancestor] [T2 swap 0 apoptosis/extrusion] ...
     *
     * where [num T2 swaps] denotes the number of T2 swaps at the present time,
     * [T2 swap 0 ancestor]  denotes the ancestor of the cell undergoing the T2
     * swap, and [T2 swap 0 apoptosis/extrusion] whether the cell was apoptotic
     * at the time of the T2 swap (literal value of 'apoptosis' if so,
     * 'extrusion' else) with index 0 in the CustomVertexBasedCellPopulation members
     * mAncestorsOfT2Swaps and mIsApoptoticsOfT2Swaps, and so on.
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the CustomVertexBasedCellPopulation to visit.
     */
    void Visit(CustomVertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(T2SwapAncestorIsApoptoticWriter)


#endif // T2SWAPANCESTORISAPOPTOTICWRITER_HPP
