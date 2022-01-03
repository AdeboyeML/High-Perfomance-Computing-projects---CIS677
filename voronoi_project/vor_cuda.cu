#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <cuda_runtime.h>
#include <cuda_runtime_api.h>

#include <iostream>
#include <iterator>
#include <vector>
#include <chrono>

__global__ void Kernel( int SizeX , int SizeY , const float2 * SiteArray , const int * Ping , int * Pong , int k , int * Mutex )
{
    //
    const int CellX = threadIdx.x + blockIdx.x * blockDim.x ;
    const int CellY = threadIdx.y + blockIdx.y * blockDim.y ;

    const int CellIdx = CellX + CellY * SizeX ;
    const int Seed = Ping[CellIdx] ;
    if ( Seed < 0 )
    {
        return ;
    }

    //
    const int2 OffsetArray[8] = { { - 1 , - 1 } ,
                                  {   0 , - 1 } ,
                                  {   1 , - 1 } ,
                                  { - 1 ,   0 } ,
                                  {   1 ,   0 } ,
                                  { - 1 ,   1 } ,
                                  {   0 ,   1 } ,
                                  {   1 ,   1 } } ;

    for ( int i = 0 ; i < 8 ; ++ i )
    {
        const int FillCellX = CellX + k * OffsetArray[i].x ;
        const int FillCellY = CellY + k * OffsetArray[i].y ; 
        if ( FillCellX >= 0 && FillCellX < SizeX && FillCellY >= 0 && FillCellY < SizeY )
        {
            //
            const int FillCellIdx = FillCellX + FillCellY * SizeX ;

            // Lock
            //
            while ( atomicCAS( Mutex , - 1 , FillCellIdx ) == FillCellIdx )
            {
            }

            const int FillSeed = Pong[FillCellIdx] ;

            if ( FillSeed < 0 )
            {
                Pong[FillCellIdx] = Seed ;
            }
            else
            {
                float2 P = make_float2( FillCellX + 0.5f , FillCellY + 0.5f ) ;

                float2 A = SiteArray[Seed] ;
                float2 PA = make_float2( A.x - P.x , A.y - P.y ) ;
                float PALength = PA.x * PA.x + PA.y * PA.y ;

                const float2 B = SiteArray[FillSeed] ;
                float2 PB = make_float2( B.x - P.x , B.y - P.y ) ;
                float PBLength = PB.x * PB.x + PB.y * PB.y ;

                if ( PALength < PBLength )
                {
                    Pong[FillCellIdx] = Seed ;
                }
            }

            // Release
            //
            atomicExch( Mutex , - 1 ) ;
        }
    }
}

int main( int Argc , char * Argv[] )
{
    -- Argc , ++ Argv ;
    if ( Argc != 3 )
    {
        printf("SOMETHING IS WRONG") ;
        return EXIT_FAILURE ;
    }

    //numSeeds - Number of Seeds
    //Size - Voronoi grid size
    int numSeeds = atoi( Argv[0] ) ;
    int Size     = atoi( Argv[1] ) ;

    //
    int NumCudaDevice = 0 ;
    cudaGetDeviceCount( & NumCudaDevice ) ;
    if ( ! NumCudaDevice )
    {
        return EXIT_FAILURE ;
    }

    //1. Generate x and y position values for the seeds in seedVec
    //2. Randomly assign seed number to some grid points(x.y) in the voronoiVec
    //3. Assign randomly generated colours to each of the seeds in randomcolourVec
    std::vector< float2 > seedVec ;
    std::vector< int >    voronoiVec( Size * Size , - 1 ) ;
    std::vector< uchar3 > randomColorVec ;
    for ( int i = 0 ; i < numSeeds ; ++ i )
    {
        float X = static_cast< float >( rand() ) / RAND_MAX * Size ;
        float Y = static_cast< float >( rand() ) / RAND_MAX * Size ;
        int CellX = static_cast< int >( floorf( X ) ) ;
        int CellY = static_cast< int >( floorf( Y ) ) ;

        seedVec.push_back( make_float2( CellX + 0.5f , CellY + 0.5f ) ) ;
        voronoiVec[CellX + CellY * Size] = i ;
        //printf("SOMETHING IS GOOD");

        randomColorVec.push_back( make_uchar3( static_cast< unsigned char >( static_cast< float >( rand() ) / RAND_MAX * 255.0f ) ,
                                               static_cast< unsigned char >( static_cast< float >( rand() ) / RAND_MAX * 255.0f ) ,
                                               static_cast< unsigned char >( static_cast< float >( rand() ) / RAND_MAX * 255.0f ) ) ) ;
    }

    //
    size_t seedSize = numSeeds * sizeof( float2 ) ;

    float2 * seedArray = NULL ;
    cudaMalloc( & seedArray , seedSize ) ;
    cudaMemcpy( seedArray , & seedVec[0] , seedSize , cudaMemcpyHostToDevice ) ;

    //BufferSize - Voronoi grid size (Size * Size)
    size_t BufferSize = Size * Size * sizeof( int ) ;

    int * Ping = NULL , * Pong = NULL ;
    cudaMalloc( & Ping , BufferSize ) , cudaMemcpy( Ping , & voronoiVec[0] , BufferSize , cudaMemcpyHostToDevice ) ;
    cudaMalloc( & Pong , BufferSize ) , cudaMemcpy( Pong , Ping , BufferSize , cudaMemcpyDeviceToDevice ) ;

    //Mutex will be used in the kernel to lock and unlock atomic operation
    int * Mutex = NULL ;
    cudaMalloc( & Mutex , sizeof( int ) ) , cudaMemset( Mutex , - 1 , sizeof( int ) ) ;

    //
    //
    cudaDeviceProp CudaDeviceProperty ;
    cudaGetDeviceProperties( & CudaDeviceProperty , 0 ) ;

    //warpsize = 32 threads
    dim3 BlockDim( CudaDeviceProperty.warpSize , CudaDeviceProperty.warpSize ) ; 
    dim3 GridDim( ( Size + BlockDim.x - 1 ) / BlockDim.x ,
                  ( Size + BlockDim.y - 1 ) / BlockDim.y ) ;

    //run JFA for logn rounds 
    auto start = std::chrono::high_resolution_clock::now();
    for ( int k = Size / 2 ; k > 0 ; k = k >> 1 )
    {
        Kernel<<< GridDim , BlockDim >>>( Size , Size , seedArray , Ping , Pong , k , Mutex ) ;
        cudaDeviceSynchronize() ;

        cudaMemcpy( Ping , Pong , BufferSize , cudaMemcpyDeviceToDevice ) ;
        std::swap( Ping , Pong ) ;
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop -start);
    printf("Execution time %ld microseconds\n", duration.count());

    cudaMemcpy( & voronoiVec[0] , Pong , BufferSize , cudaMemcpyDeviceToHost ) ;

    //
    cudaFree( seedArray ) ;
    cudaFree( Ping ) ;
    cudaFree( Pong ) ;
    cudaFree( Mutex ) ;

    //
    //
    FILE * Output = fopen( Argv[2], "wb" ) ;
    fprintf( Output , "P6\n%d %d\n255\n" , Size , Size ) ;

    std::vector< uchar3 > Pixels( Size * Size ) ;
    for ( int y = 0 ; y < Size ; ++ y )
    {
        for ( int x = 0 ; x < Size ; ++ x )
        {
            const int Seed = voronoiVec[x + y * Size] ;
            if ( Seed != - 1 )
            {
                Pixels[x + y * Size] = randomColorVec[Seed] ;
            }
        }
    }

    for( std::vector< float2 >::const_iterator itr = seedVec.begin() ; itr != seedVec.end() ; ++ itr )
    {
        const int x = static_cast< int >( floorf( itr->x ) ) ;
        const int y = static_cast< int >( floorf( itr->y ) ) ;
        Pixels[x + y * Size] = make_uchar3( 255 , 0 , 0 ) ;
    }

    fwrite( & Pixels[0].x , 3 , Pixels.size() , Output ) ;
    fclose( Output ) ;

    return EXIT_SUCCESS ;
}