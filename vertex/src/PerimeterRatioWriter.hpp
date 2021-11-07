#ifndef PERIMETERRATIOWRITER_HPP
#define PERIMETERRATIOWRITER_HPP

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


/**
 * A class written using the visitor pattern for writing the proportion of the
 * perimeter that is bordered by a cell of the same ancestor, weighted by cell
 * areas.
 *
 * The output file is called perimeterratio.dat by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class PerimeterRatioWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
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
    PerimeterRatioWriter();

    /**
     * Visit the population and write the perimeter ratio data.
     *
     * Outputs a line of tab-separated values of the form:
     * ...[cell_ancestor_index] [perimeter_ratio] ...
     *
     * Here the indexing of CellAncestor is as given by Cell::GetAncestor().
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation) override;

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is not defined for use
     * with a CaBasedCellPopulation.
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation) override;

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is not defined for use
     * with a MeshBasedCellPopulation.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation) override;

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is not defined for use
     * with a NodeBasedCellPopulation.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation) override;

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is not defined for use
     * with a PottsBasedCellPopulation.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation) override;
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(PerimeterRatioWriter)

#endif // PERIMETERRATIOWRITER_HPP
