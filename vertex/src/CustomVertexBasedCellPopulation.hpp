#ifndef CUSTOMVERTEXBASEDCELLPOPULATION_HPP
#define CUSTOMVERTEXBASEDCELLPOPULATION_HPP

#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
class CustomVertexBasedCellPopulation : public VertexBasedCellPopulation<DIM>
{
private:

    /**
     * Ancestors of T2 swaps (the centre of the removed triangle), stored so
     * they can be accessed and output by the cell killer and
     * population writer classes. The locations are stored until they are
     * cleared by ClearAncestorsAndIsApoptoticsOfT2Swaps();
     */
    std::vector<unsigned> mAncestorsOfT2Swaps;

    /**
     * Flag indicating whether the cells that have undergone T2 swaps were
     * apoptotic, stored so they can be accessed and output by the cell killer
     * and population writer classes. The locations are stored until they are
     * cleared by ClearAncestorsAndIsApoptoticsOfT2Swaps();
     */
    std::vector<bool> mIsApoptoticsOfT2Swaps;

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<VertexBasedCellPopulation<DIM> >(*this);
    }

public:

    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely one CellPtr for each VertexElement in
     * the mesh.
     *
     * @param rMesh reference to a
     * @param rCells reference to a vector of CellPtrs
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     * @param validate whether to validate the cell population when it is created (defaults to true)
     * @param locationIndices an optional vector of location indices that correspond to real cells
     */
    CustomVertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh,
                              std::vector<CellPtr>& rCells,
                              bool deleteMesh=false,
                              bool validate=true,
                              const std::vector<unsigned> locationIndices=std::vector<unsigned>());

    /**
     * Constructor for use by boost serialization ONLY!
     *
     * @param rMesh a vertex mesh.
     */
    CustomVertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~CustomVertexBasedCellPopulation() override;

    /**
     * Return ancestors of T2 swaps since the last sampling time step.
     *
     * @return mAncestorsOfT2Swaps
     */
    std::vector<unsigned> GetAncestorsOfT2Swaps();

    /**
     * Return vector of booleans indicating whether the cells undergoing T2
     * swaps were apoptotic at the time of the T2 swap since the last sampling
     * time step.
     *
     * @return mIsApoptoticsOfT2Swaps
     */
    std::vector<bool> GetIsApoptoticsOfT2Swaps();

    /**
     * Add ancestor of T2 Swap.
     *
     * @param ancestorOfT2Swap  Ancestor of the T2 swap
     */
    void AddAncestorOfT2Swap(unsigned ancestorOfT2Swap);

    /**
     * Add boolean indicating whether cell undergoing T2 Swap was apoptotic at
     * the time of the T2 swap
     *
     * @param isApoptoticOfT2Swap  whether cell was apopototic
     */
    void AddIsApoptoticOfT2Swap(bool isApoptoticOfT2Swap);

    /**
     * Clear the ancestors and is apoptotic flags of cells undergoing T2 Swaps.
     */
    void ClearAncestorsAndIsApoptoticsOfT2Swaps();

    /**
     * Overridden SimulationSetupHook() method.
     *
     * Hook method to add a CustomT2SwapCellKiller to a simulation object, which is always
     * required in the case of a VertexBasedCellPopulation. This functionality avoids
     * the need for static or dynamic casts to specific cell population types within
     * simulation methods.
     *
     * Note: In order to inhibit T2 swaps, the user needs to set the threshold for T2
     * swaps in the MutableVertexMesh object mrMesh to 0, using the SetT2Threshold() method.
     *
     * @param pSimulation pointer to a cell-based simulation object
     */
    virtual void SimulationSetupHook(AbstractCellBasedSimulation<DIM, DIM>* pSimulation) override;

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CustomVertexBasedCellPopulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a CustomVertexBasedCellPopulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const CustomVertexBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const MutableVertexMesh<DIM,DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a CustomVertexBasedCellPopulation.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CustomVertexBasedCellPopulation<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    MutableVertexMesh<DIM,DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)CustomVertexBasedCellPopulation<DIM>(*p_mesh);
}
}
} // namespace ...
#endif // CUSTOMVERTEXBASEDCELLPOPULATION_HPP
