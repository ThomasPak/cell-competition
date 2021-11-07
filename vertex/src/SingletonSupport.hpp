#ifndef SINGLETONSUPPORT_HPP_
#define SINGLETONSUPPORT_HPP

/** Static methods for setting up and destroying singletons.
 *
 * The singleton functions are taken from AbstractCellBasedTestSuite::setUp()
 * and AbstractCellBasedTestSuite::tearDown().
 */

class SingletonSupport
{
    public:

        /** Setup singletons.
         *
         * Seed random number generator with seed.  If seed is 0,
         * use time as seed.
         *
         * @param seed  seed for random number generator.
         */
        static void SetupSingletons(unsigned seed = 0);

        /** Destroy singletons. */
        static void DestroySingletons();
};

#endif // SINGLETONSUPPORT_HPP_
