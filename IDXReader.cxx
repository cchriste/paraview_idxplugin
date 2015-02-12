#include "IDXReader.h"
 
#include <vtkShortArray.h>
#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkStructuredGrid.h>
#include <vtkImageData.h>
#include <vtkRectilinearGrid.h>

#include "vtkObjectFactory.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkInformation.h"
#include "vtkDataObject.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCharArray.h"

#ifdef AM   /* already included in IDXReader.h */
#include <string>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>

#include <visuscpp/idx/visus_idx.h>
#endif

using std::string;
VISUS_USE_NAMESPACE

UniquePtr<IdxDataset> dataset;
//UniquePtr<Dataset> dataset;
//ScopedPtr<Dataset> dataset;
UniquePtr<Application> app;
//ScopedPtr<Application> app;
UniquePtr<Access> accessdata;
//ScopedPtr<Access> accessdata;

vtkStandardNewMacro(IDXReader);
 
/********************************************************************************/
IDXReader::IDXReader()
{
	string fieldname = "DATA";

	std::cout << "IDXReader: IDXReader() ..." << endl;
  this->FileName = NULL;
  this->PercentToRemove = 4;

	this->FieldName = strdup(fieldname.c_str());
	//AM std::cout << "  this->FieldName=" << this->FieldName << endl;

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

/********************************************************************************/
void IDXReader::SetPercentToRemove(double percent)
{
	std::cout << "IDXReader: SetPercentToRemove ..." << endl;
	this->PercentToRemove = percent;

	this->Modified();
}

/********************************************************************************/
//AM void IDXReader::SetFieldName(char* field_name)
//AM {
//AM 	std::cout << "IDXReader: SetFieldName ..." << endl;
//AM 
//AM 	//std::cout << "  field_name=" << field_name << endl;
//AM 	this->FieldName = field_name;
//AM 	std::cout << "  this->FieldName=" << this->FieldName << endl;
//AM 
//AM 	this->Modified();
//AM }

/********************************************************************************/
void IDXReader::SetTransformMode(int mode)
{
	std::cout << "IDXReader: SetTransformMode ..." << endl;
	this->TransformMode = mode;
}


/********************************************************************************/
int IDXReader::RequestInformation(
    vtkInformation * vtkNotUsed(request),
    vtkInformationVector ** vtkNotUsed(inputVector),
    vtkInformationVector *outputVector)
{
   int testingnew=0;
    
   std::cout << "IDXReader: RequestInformation ..." << endl;
   vtkInformation* outInfo = outputVector->GetInformationObject(0);

   /*  
   int extent[6];
 
   extent[0] = 0;
   extent[1] = 127;
   extent[2] = 0;
   extent[3] = 127;
   extent[4] = 0;
   extent[5] = 127;
 
   outInfo->Set
     (vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent,6);
   */

   std::cerr << "IDXReader: requesting information return" << endl;
   return 1;
}


/********************************************************************************/
int IDXReader::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **vtkNotUsed(inputVector),
  vtkInformationVector *outputVector)
{
 
	//std::cout << "@@@@@  IDXReader: RequestData ..." << endl;
  
  /*
  int extents[3];
  int bounds[3];
  */
  int slice; 
  int sub_div[3];
  int proc[3];
  int *count_local;
  int *r_count_local;
  int *offset_local;
  int *r_offset_local;
	//vtkRectilinearGrid* img;
  vtkRectilinearGrid *img = vtkRectilinearGrid::New();

	//SharedPtr<Dataset> dataset;
	//UniquePtr<Dataset> dataset;
	//ScopedPtr<Dataset> dataset;
	//UniquePtr<Application> app;
	//ScopedPtr<Application> app;
	//UniquePtr<Access> accessdata;
	//ScopedPtr<Access> accessdata;

	static int bounds[3];
	static double extents[6];

  string name(this->FileName);
  int resolution = (int)this->PercentToRemove;
  int rank, nprocs;
  int domain, nblocks;
  int url_open = 0;

  // get the info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
  /*
  int subext[6];
  outInfo->Get
      (vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(),subext);
  
  std::cerr << "PER PROCESS=("
       << subext[0] << "," << subext[1] << "),("
       << subext[2] << "," << subext[3] << "),("
       << subext[4] << "," << subext[5] << ")" << endl; 
  */

  rank = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());

  nprocs =
    outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  
  nblocks = nprocs;
  output->SetNumberOfBlocks(1);
  vtkRectilinearGrid* rgrid = vtkRectilinearGrid::New();
  rgrid->Initialize();
  output->SetBlock(0,rgrid);

  //ScopedPtr<Application> app(new Application);
  app.reset(new Application);
  app->setCommandLine(0, NULL);
  ENABLE_VISUS_IDX();
 	//std::cout << "IDXReader: RESOLUTION, PercentToRemove=" << this->PercentToRemove << endl;
  if(url_open == 1)
  	name = "http://atlantis.sci.utah.edu/mod_visus?action=readdataset&dataset=2kbit1";
  else
		name = (this->FileName);
  //std::cout << "IDXReader: name to open is " << name << endl;
	
#ifdef AM
  dataset.reset(Dataset::open(name));
#else
  //dataset.reset(Dataset::loadDataset(name));
  dataset.reset(dynamic_cast<IdxDataset*>(Dataset::loadDataset(name)));
#endif
  VisusReleaseAssert(dataset);
    
  //Redundant Code
  NdBox world_box = dataset->logic_box;
  //std::cout << "IDXReader: after NdBox world_box" << endl;

  accessdata.reset(dataset->createAccess());

  bounds[0] = world_box.to.x - world_box.from.x + 1;
  bounds[1] = world_box.to.y - world_box.from.y + 1;
  bounds[2] = world_box.to.z - world_box.from.z + 1;
  cerr << "IDXReader: dataset->logic_box=("
       << world_box.from.x << "," << world_box.to.x << "),("
       << world_box.from.y << "," << world_box.to.y << "),("
       << world_box.from.z << "," << world_box.to.z << ")" << endl;
  cerr << "IDXReader: dataset->maxh=" << dataset->max_resolution << endl;
    
  // Convert the logical extents into physical extents.
  extents[0] = world_box.from.x;
  extents[1] = world_box.to.x + 1;

  extents[2] = world_box.from.y;
  extents[3] = world_box.to.y + 1;

  extents[4] = world_box.from.z;
  extents[5] = world_box.to.z + 1;

  NdBox slice_box = dataset->logic_box;

  //cerr << "IDXReader: GetMesh:resolution=" << resolution << endl;

  //
  // Determine the reduced grid size taking into account the multi-res
  // reduction.
  //
  int res = resolution;
  int resReduction = 1 << res;
  //cerr << "IDXReader: resReduction=" << resReduction << endl;
  int reducedBounds[3];

  domain = rank;

  reducedBounds[0] = bounds[0] % resReduction == 0 ?
      bounds[0] / resReduction : bounds[0] / resReduction + 1;
  reducedBounds[1] = bounds[1] % resReduction == 0 ?
      bounds[1] / resReduction : bounds[1] / resReduction + 1;
  reducedBounds[2] = bounds[2] % resReduction == 0 ?
      bounds[2] / resReduction : bounds[2] / resReduction + 1;

  int zWidth = reducedBounds[2] % nblocks == 0 ?
        reducedBounds[2] / nblocks : reducedBounds[2] / nblocks + 1;
	if(reducedBounds[2] % nblocks == 0)
	{
		slice_box.from[2] = zWidth * resReduction * domain;
		slice_box.to[2]   = slice_box.from[2] + zWidth * resReduction - 1;
	}
	else
	{
		if(domain == nblocks - 1)
		{
			slice_box.from[2] = zWidth * resReduction * domain;
			slice_box.to[2]   = slice_box.to[2];
		}
		else
		{       
			slice_box.from[2] = zWidth * resReduction * domain;
			slice_box.to[2]   = slice_box.from[2] + zWidth * resReduction - 1;
		}
  }

    //
    // Create the mesh.
    //
    //vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
		/* */
/* Block 1 commented out */
/*  */
#ifdef block1_commented_out
#endif
    int dims[3];
    float *arrayX;
    float *arrayZ;
    float *arrayY;
    vtkFloatArray *coordsX;
    vtkFloatArray *coordsY;
    vtkFloatArray *coordsZ;
    
    dims[0] = (slice_box.to[0] - slice_box.from[0] ) / resReduction + 1;
    dims[1] = (slice_box.to[1] - slice_box.from[1] ) / resReduction + 1;
    dims[2] = (slice_box.to[2] - slice_box.from[2] ) / resReduction + 1;
    rgrid->SetDimensions(dims[0], dims[1], dims[2]);

    coordsX = vtkFloatArray::New();
    coordsX->SetNumberOfTuples(dims[0]);
    arrayX = (float *) coordsX->GetVoidPointer(0);
    for (int i = 0; i < dims[0]; i++)
        arrayX[i] = i * resReduction;
    arrayX[dims[0]-1] = slice_box.to[0] + 1.;
    rgrid->SetXCoordinates(coordsX);

    coordsY = vtkFloatArray::New();
    coordsY->SetNumberOfTuples(dims[1]);
    arrayY = (float *) coordsY->GetVoidPointer(0);
    for (int i = 0; i < dims[1]; i++)
        arrayY[i] = i * resReduction;
    arrayY[dims[1]-1] = slice_box.to[1] + 1.;
    rgrid->SetYCoordinates(coordsY);

    coordsZ = vtkFloatArray::New();
    coordsZ->SetNumberOfTuples(dims[2]);
    arrayZ = (float *) coordsZ->GetVoidPointer(0);
    for (int i = 0; i < dims[2]; i++)
        arrayZ[i] = slice_box.from[2] + i * resReduction;
    arrayZ[dims[2]-1] = slice_box.to[2] + 1.;
    rgrid->SetZCoordinates(coordsZ);
    

    Field field; 
    // -AM 20140327  // field = dataset->getFieldByName("var");
		//
		for (std::vector<Field>::iterator it=dataset->idxfile.fields.begin();it!=dataset->idxfile.fields.end();it++)
		{
			field = (*it);
			break;
		}

		if (field.valid())
		{
		  //std::cout << "IDXReader: RequestData: valid field: field.name=" << field.name << endl;
		  this->FieldName = strdup(field.name.c_str());
		}
		else
		{
			std::cout << "IDXReader: RequestData: *********  field is NOT valid " << endl;
		}

    
    int maxh_read = dataset->max_resolution - res * 3;
    
#ifdef AM
    //BoxQuery specs(timestat/0, dataset->getFieldByName("var"),Position::LogicSpace,slice_box);
    BoxQuery specs(0, dataset->getFieldByName("var"),Position::LogicSpace,slice_box);
#else
    //BoxQuerySpec specs(0, dataset->getFieldByName("var"), Position::LogicSpace, slice_box);
    BoxQuerySpec specs(0, field, Position::LogicSpace, slice_box);
#endif

    specs.setResolutionRange(0, maxh_read, dataset->max_resolution);
    
#ifdef AM
    //BoxQueryPtr query = dataset->createBoxQuery(Position::LogicSpace,
    //                    slice_box, dataset->getFieldByName(name), timestate, 0,
    //                    maxh_read, dataset->max_resolution);
    ScopedPtr<BoxQueryJob> query(dataset->createJob(specs));
#else
    UniquePtr<BoxQuery> query(dataset->createBoxQuery(specs));
#endif
		//std::cout << "IDXReader: after uniqueptr BoxQuery query" << endl;

		int ntuples;
#ifdef AM
		ntuples = query->nsamples.innerProduct();
		VisusReleaseAssert(query);
#else
		ntuples = query->nsamples.innerProduct();
		//std::cout << "IDXReader: ntuples=" << ntuples << endl;
		if (specs.valid() )
		{
			VisusReleaseAssert(query);
			//std::cout << "IDXReader: valid, after assert " << endl;
		}
#endif


/* Block 2 commented out */
/*  */
#ifdef block2_commented_out
#endif

		//std::cout << "IDXReader: before DType if statement" << endl;
    if (field.dtype == DType::UINT8)
    {
        std::cout << "DType=UINT8"  << endl;

        vtkUnsignedCharArray *rv1 = vtkUnsignedCharArray::New();
        
        //rv1->SetNumberOfComponents(1);
        //rv1->SetNumberOfTuples(ntuples);
        
        //VisusReleaseAssert(rv1->GetPointer(0) != NULL);
        //query->buffer=SharedPtr<Array>(new Array((unsigned char*) rv1->GetPointer(0), query->nsamples, DType::FLOAT64));
        //std::cout << "[1]" << query->buffer->c_size() << endl;
        //std::cout << "[2]" << sizeof (double) * ntuples << endl;
	
        VisusReleaseAssert(query->allocateBuffers());
#ifdef AM
        dataset->execute(query.get(), accessdata.get(),'r');
#else
        dataset->executeBoxQuery(query.get(), accessdata.get(),'r');
#endif

        //std::cout << "query buffer c_size=" << query->buffer->c_size() << endl;
        //std::cout << "ntuples=" << sizeof (unsigned char) * ntuples << endl;

        VisusReleaseAssert(query->ok());
        //VisusReleaseAssert(query->buffer->c_size()==sizeof(double)*rcount[0]*rcount[1]*rcount[2]);
        
        VisusReleaseAssert(query->buffer->c_size() == sizeof (unsigned char) * ntuples);
        
        //dataset->execute(query.get(), accessdata.get(),'r');
				//VisusReleaseAssert(query->ok());
        
        unsigned char* Src=(unsigned char*)query->buffer->c_ptr();
#ifdef AM
        for (int I=0;I<rcount[0]*rcount[1]*rcount[2] ntuples;I++)
#else
        for (int I=0;I<ntuples;I++)
#endif
        {
           //if(I < 128*128*128 && I%128 == 0)
					 // std::cerr<< I << " : " << Src[I] << endl;

           rv1->InsertNextValue(Src[I]);
        }

        rv1->SetName("data");
        rgrid->GetPointData()->SetScalars(rv1);

				//adding data to rectilineargrid
	

    }
    else if (field.dtype == DType::FLOAT64)
    {
        //std::cout << "DType=FLOAT64"  << endl;

        vtkDoubleArray *rv1 = vtkDoubleArray::New();
        //rv1->SetNumberOfComponents(1);
        //rv1->SetNumberOfTuples(ntuples);
        
        //VisusReleaseAssert(rv1->GetPointer(0) != NULL);
        //query->buffer=SharedPtr<Array>(new Array((unsigned char*) rv1->GetPointer(0), query->nsamples, DType::FLOAT64));
        //std::cout << "[1]" << query->buffer->c_size() << endl;
        //std::cout << "[2]" << sizeof (double) * ntuples << endl;
	
        VisusReleaseAssert(query->allocateBuffers());
#ifdef AM
        dataset->execute(query.get(), accessdata.get(),'r');
#else
        dataset->executeBoxQuery(query.get(), accessdata.get(),'r');
#endif
        //std::cout << "query buffer c_size=" << query->buffer->c_size() << endl;
        //std::cout << "ntuples=" << sizeof (unsigned char) * ntuples << endl;
        VisusReleaseAssert(query->ok());

        //VisusReleaseAssert(query->buffer->c_size()==sizeof(double)*rcount[0]*rcount[1]*rcount[2]);
        
        VisusReleaseAssert(query->buffer->c_size() == sizeof (double) * ntuples);
        
        //dataset->execute(query.get(), accessdata.get(),'r');
				//VisusReleaseAssert(query->ok());
        
        double* Src=(double*)query->buffer->c_ptr();
#ifdef AM
        for (int I=0;I<rcount[0]*rcount[1]*rcount[2] ntuples;I++)
#else
        for (int I=0;I<ntuples;I++)
#endif
        {
           //if(I < 128*128*128 && I%128 == 0)
					 // std::cerr<< I << " : " << Src[I] << endl;

           rv1->InsertNextValue(Src[I]);
        }

        rv1->SetName("data");
        rgrid->GetPointData()->SetScalars(rv1);

        //return rv1;
    }
    else if (field.dtype == DType::FLOAT32)
    {
        //std::cout << "DType=float 32 "  << endl;

        vtkFloatArray *rv1 = vtkFloatArray::New();
        //std::cout << "after vtkfloatarray "  << endl;
        VisusReleaseAssert(query->allocateBuffers());
        //std::cout << "after releaseAssert "  << endl;
#ifdef AM
        dataset->execute(query.get(), accessdata.get(),'r');
#else
        dataset->executeBoxQuery(query.get(), accessdata.get(),'r');
        //std::cout << "after execute "  << endl;
#endif

        //std::cout << "query buffer c_size=" << query->buffer->c_size() << endl;
        //std::cout << "ntuples=" << sizeof (unsigned char) * ntuples << endl;
        VisusReleaseAssert(query->ok());
        //std::cout << "after releaseassert ok"  << endl;

        //VisusReleaseAssert(query->buffer->c_size()==sizeof(double)*rcount[0]*rcount[1]*rcount[2]);
        
        VisusReleaseAssert(query->buffer->c_size() == sizeof (float) * ntuples);
        //std::cout << "after releaseassert buffer"  << endl;
        
        //dataset->execute(query.get(), accessdata.get(),'r');
				//VisusReleaseAssert(query->ok());
        
        float* Src=(float*)query->buffer->c_ptr();
        //std::cout << "Src="  << *Src << endl;

#ifdef AM
        for (int I=0;I<rcount[0]*rcount[1]*rcount[2] ntuples;I++)
#else
        for (int I=0;I<ntuples;I++)
#endif
        {
           //if(I < 128*128*128 && I%128 == 0)
					 // std::cerr<< I << " : " << Src[I] << endl;

           rv1->InsertNextValue(Src[I]);
        }

        rv1->SetName("data");
        rgrid->GetPointData()->SetScalars(rv1);

	
        //return rv;
    }
    else if (field.dtype == DType::UINT16)
    {
        //std::cout << "DTYPE=UINT16 "  << endl;
				//std::cout << "DType is:" << field.dtype << endl;

        vtkShortArray *rv1 = vtkShortArray::New();
        rv1->SetNumberOfComponents(1);
        rv1->SetNumberOfTuples(ntuples);

				VisusReleaseAssert(rv1->GetPointer(0) != NULL);
        query->buffer=SharedPtr<Array>(new Array((unsigned char*) rv1->GetPointer(0), query->nsamples, DType::UINT16));
				VisusReleaseAssert(query->buffer->c_size() == sizeof (short) * ntuples);
	
#ifdef AM
        dataset->execute(query.get(), accessdata.get(),'r');
#else
        dataset->executeBoxQuery(query.get(), accessdata.get(),'r');
#endif
				VisusReleaseAssert(query->ok());
        
        //return rv1;
    }
		else 
		{
			std::cout << "else on DType is ???" << endl;
			//std::cout << "DType is:" << field.dtype << endl;
		}

/* Block 3 commented out */
/*  */

#ifdef block3_commented_out
#endif

#ifdef block3_tmp
  app.reset(new Application);
  app->setCommandLine(0, NULL);
  ENABLE_VISUS_IDX();  

  if(url_open == 1)
      name = "http://atlantis.sci.utah.edu/mod_visus?action=readdataset&dataset=2kbit1";
  else
      name = (this->FileName);
  //std::cout << "IDXReader: name in block3 to open is " << name << endl;
 
  //cerr << "[1]" << endl;
 
#ifdef AM
  dataset.reset(Dataset::open(name));
#else
  //dataset.reset(Dataset::loadDataset(name));
  dataset.reset(dynamic_cast<IdxDataset*>(Dataset::loadDataset(name)));
#endif
  VisusReleaseAssert(dataset);

  //cerr << "[2]" << endl;

  //Redundant Code
#ifdef AM
  NdBox world_box = dataset->logic_box;
#endif
    
  //cerr << "[3]" << endl;

  accessdata.reset(dataset->createAccess());

  //cerr << "[4]" << endl;    

  bounds[0] = world_box.to.x - world_box.from.x + 1;
  bounds[1] = world_box.to.y - world_box.from.y + 1;
  bounds[2] = world_box.to.z - world_box.from.z + 1;
  cerr << "dataset->logic_box=("
       << world_box.from.x << "," << world_box.to.x << "),("
       << world_box.from.y << "," << world_box.to.y << "),("
       << world_box.from.z << "," << world_box.to.z << ")" << endl;
  cerr << "dataset->maxh=" << dataset->max_resolution << endl;
    
  // Convert the logical extents into physical extents.
  extents[0] = world_box.from.x;
  extents[1] = world_box.to.x + 1;

  extents[2] = world_box.from.y;
  extents[3] = world_box.to.y + 1;

  extents[4] = world_box.from.z;
  extents[5] = world_box.to.z + 1;
#endif        /*  #ifdef block3_tmp  */
	
#ifdef AM
  int res = 3;
  int resReduction = 1 << res;
  int reducedBounds[3];
#else
  res = 3;
  resReduction = 1 << res;
#endif


  reducedBounds[0] = bounds[0] % resReduction == 0 ?
      bounds[0] / resReduction : bounds[0] / resReduction + 1;
  reducedBounds[1] = bounds[1] % resReduction == 0 ?
      bounds[1] / resReduction : bounds[1] / resReduction + 1;
  reducedBounds[2] = bounds[2] % resReduction == 0 ?
      bounds[2] / resReduction : bounds[2] / resReduction + 1;


  //extents[0] = 128;
  //extents[1] = 128;
  //extents[2] = 128;

  count_local = (int*)malloc(sizeof(int) * 3);
  offset_local = (int*)malloc(sizeof(int) * 3);

  r_count_local = (int*)malloc(sizeof(int) * 3);
  r_offset_local = (int*)malloc(sizeof(int) * 3);

  proc[0] = 2;
  proc[1] = 2;
  proc[2] = 2;

  r_count_local[0] = reducedBounds[0]/proc[0];
  r_count_local[1] = reducedBounds[1]/proc[1];
  r_count_local[2] = reducedBounds[2]/proc[2];

  //std::cout << "IDXReader: reducedBounds[0]=" << reducedBounds[0] << endl;
  //std::cout << "IDXReader: r_count_local[0]=" << r_count_local[0] << endl;
  sub_div[0] = (reducedBounds[0] / r_count_local[0]);

  //std::cout << "IDXReader: reducedBounds[1]=" << reducedBounds[1] << endl;
  //std::cout << "IDXReader: r_count_local[1]=" << r_count_local[1] << endl;
  sub_div[1] = (reducedBounds[1] / r_count_local[1]);

  std::cout << "IDXReader: reducedBounds[2]=" << reducedBounds[2] << endl;
  std::cout << "IDXReader: r_count_local[2]=" << r_count_local[2] << endl;
	if ( r_count_local[2] != 0 )
		sub_div[2] = (reducedBounds[2] / r_count_local[2]);

  r_offset_local[2] = (rank / (sub_div[0] * sub_div[1])) * r_count_local[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  std::cout << "IDXReader: after slice" << endl;
  r_offset_local[1] = (slice / sub_div[0]) * r_count_local[1];
  std::cout << "IDXReader: after r_offset_local" << endl;
  r_offset_local[0] = (slice % sub_div[0]) * r_count_local[0];

  cerr << rank << endl;
  cerr << "[LOCAL] dataset->logic_box=("
       << r_offset_local[0] << "," << r_count_local[0] << "),("
       << r_offset_local[1] << "," << r_count_local[1] << "),("
       << r_offset_local[2] << "," << r_count_local[2] << ")" << endl;


  count_local[0] = bounds[0]/proc[0];
  count_local[1] = bounds[1]/proc[1];
  count_local[2] = bounds[2]/proc[2];

  sub_div[0] = (bounds[0] / count_local[0]);
  sub_div[1] = (bounds[1] / count_local[1]);
  sub_div[2] = (bounds[2] / count_local[2]);
  offset_local[2] = (rank / (sub_div[0] * sub_div[1])) * count_local[2];
  slice = rank % (sub_div[0] * sub_div[1]);
  offset_local[1] = (slice / sub_div[0]) * count_local[1];
  offset_local[0] = (slice % sub_div[0]) * count_local[0];

  cerr << rank << endl;
  cerr << "[LOCAL] dataset->logic_box=("
       << offset_local[0] << "," << count_local[0] << "),("
       << offset_local[1] << "," << count_local[1] << "),("
       << offset_local[2] << "," << count_local[2] << ")" << endl;

  
  // Here is where you would read the data from the file. 
  ReadFile(this->FieldName, r_offset_local, r_count_local, offset_local, count_local, img);

	std::cout << "@@@@@   IDXReader: RequestData return " << endl;

  return 1;
}
 

/********************************************************************************/
void IDXReader::PrintSelf(ostream& os, vtkIndent indent)
{
	std::cout << "IDXReader: PrintSelf ..." << endl;
  this->Superclass::PrintSelf(os,indent);
 	printf("IDXReader: PrintSelf FileName: %s\n", this->FileName);
  os << indent << "File Name: "
      << (this->FileName ? this->FileName : "(none)") << "\n";
}


/********************************************************************************/
int IDXReader::ReadFile(char* field_name, int* roffset, int* rcount, int* offset, int* count, vtkRectilinearGrid* output)
{
  //int url_open = 0;
  //string name(this->FileName);  
  //static int bounds[3];

  //ScopedPtr<Dataset> dataset;
	//UniquePtr<Dataset> dataset;
  //ScopedPtr<Access> accessdata;
	//UniquePtr<Access> accessdata;
	//not needed//UniquePtr<Application> app;

  static double extents[6];
  int domain = 0;
  int nblocks = 1;
  int dims[3];

	std::cout << "#### IDXReader: ReadFile ..." << endl;


  //ViSUS Initialization
  //ScopedPtr<Application> app(new Application);
  //app->setCommandLine(0, NULL);
  //ENABLE_VISUS_IDX();

  //if(url_open == 1) name = "http://atlantis.sci.utah.edu/mod_visus?action=readdataset&dataset=2kbit1";
  //else name = (this->FileName);
 
  
  //dataset.reset(Dataset::open(name));
  //VisusReleaseAssert(dataset);
  //NdBox world_box = dataset->logic_box;
  //accessdata.reset(dataset->createAccess());
    
  //bounds[0] = world_box.to.x - world_box.from.x + 1;
  //bounds[1] = world_box.to.y - world_box.from.y + 1;
  //bounds[2] = world_box.to.z - world_box.from.z + 1;
  //cerr << "dataset->logic_box=("
  //     << world_box.from.x << "," << world_box.to.x << "),("
  //     << world_box.from.y << "," << world_box.to.y << "),("
  //     << world_box.from.z << "," << world_box.to.z << ")" << endl;
  //cerr << "dataset->maxh=" << dataset->max_resolution << endl;
    
  // Convert the logical extents into physical extents.
  //extents[0] = world_box.from.x;
  //extents[1] = world_box.to.x + 1;

  //extents[2] = world_box.from.y;
  //extents[3] = world_box.to.y + 1;

  //extents[4] = world_box.from.z;
  //extents[5] = world_box.to.z + 1;
  
 

  NdBox slice_box;// = dataset->logic_box;
  slice_box.from[0]=offset[0];
  slice_box.to  [0]=offset[0] + count[0] - 1;

  slice_box.from[1]=offset[1];
  slice_box.to  [1]=offset[1] + count[1] - 1;

  slice_box.from[2]=offset[2];
  slice_box.to  [2]=offset[2] + count[2] - 1;


  //dims[0] = (slice_box.to[0] - slice_box.from[0] + 1);
  //dims[1] = (slice_box.to[1] - slice_box.from[1] + 1);
  //dims[2] = (slice_box.to[2] - slice_box.from[2] + 1);
  //output->SetDimensions(dims[0], dims[1], dims[2]);
  output->SetDimensions(rcount[0], rcount[1], rcount[2]);

  //std::cerr << "IDXReader: D0 " << dims[0] << endl;
  //std::cerr << "IDXReader: D1 " << dims[1] << endl;
  //std::cerr << "IDXReader: D2 " << dims[2] << endl;

  //vtk coordinate setup
  vtkFloatArray *coordsX = vtkFloatArray::New();
  for (float i = roffset[0]; i < roffset[0] + rcount[0]; i++)
		coordsX->InsertNextValue(i); 
  output->SetXCoordinates(coordsX);

  vtkFloatArray *coordsY = vtkFloatArray::New();
  for (float j = roffset[1]; j < roffset[1] + rcount[1]; j++)
    coordsY->InsertNextValue(j); 
  output->SetYCoordinates(coordsY);

  vtkFloatArray *coordsZ = vtkFloatArray::New();  
  for (float k = roffset[2]; k < roffset[2] + rcount[2]; k++)
    coordsZ->InsertNextValue(k); 
  output->SetZCoordinates(coordsZ);

  //ViSUS data read
  
  Field field; 
	// -AM 20140327  // field = dataset->getFieldByName("var");
	field = dataset->getFieldByName(field_name);
	if (field.valid())
	{
		std::cout << "ReadFile: valid field: field.name=" << field.name << endl;
	}
	else
	{
		std::cout << "  *********  field is NOT valid " << endl;
	}

  int maxh_read = dataset->max_resolution - /*res*/3 * 3;

#ifdef AM
  BoxQuery specs(0, dataset->getFieldByName("var"), Position::LogicSpace, slice_box);
#else
  //BoxQuerySpec specs(0, dataset->getFieldByName("var"), Position::LogicSpace, slice_box);
  BoxQuerySpec specs(0, field, Position::LogicSpace, slice_box);
#endif
  specs.setResolutionRange(0, maxh_read, dataset->max_resolution);
    
 
#ifdef AM
  ScopedPtr<BoxQueryJob> query(dataset->createJob(specs));
#else
  UniquePtr<BoxQuery> query(dataset->createBoxQuery(specs));
#endif
  int ntuples = query->nsamples.innerProduct();
  VisusReleaseAssert(query);

  VisusReleaseAssert(query->allocateBuffers());
#ifdef AM
  dataset->execute(query.get(), accessdata.get(), 'r');
#else
  dataset->executeBoxQuery(query.get(), accessdata.get(), 'r');
#endif
  VisusReleaseAssert(query->ok());
  std::cerr << "IDXReader: [1] query->buffer->c_size=" << query->buffer->c_size() << endl;
  std::cerr << "IDXReader: [2] rcount0 * rcount1 * rcount2=" << sizeof(double)*rcount[0]*rcount[1]*rcount[2] << endl;
  std::cerr << "IDXReader: [3] ntuples=" << ntuples  << endl;
  // -AM//VisusReleaseAssert(query->buffer->c_size()==sizeof(double)*rcount[0]*rcount[1]*rcount[2]);
 

  //adding data to rectilineargrid
  vtkDoubleArray* newArray = newArray = vtkDoubleArray::New();
  double* Src=(double*)query->buffer->c_ptr();


	// hard code values for data  -- temp
	//srand(time(NULL));
	//int anum;
  //for (int I=0;I<rcount[0]*rcount[1]*rcount[2];I++)
	//for (int I=0;I<ntuples;I++)
  //{
	//	anum = rand() % 300;
	//	Src[I] =  anum / 100;
  //
  //		if ( I < 10 )
	//		std::cerr << "IDXReader: Src[I]=" << Src[I] << endl;
  //
	//	newArray->InsertNextValue(Src[I]);
	//}


  newArray->SetName("DATA");
  output->GetPointData()->SetScalars(newArray);
  
	std::cout << "#### IDXReader: ReadFile return " << endl;

  return 1;//28*128*128;//ntuples;

}
