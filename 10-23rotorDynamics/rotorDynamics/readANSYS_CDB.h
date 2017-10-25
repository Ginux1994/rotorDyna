/*
*	main.c
*/
#include "stdafx.h"
#include  <stdio.h>
#include  <math.h>
#include <limits.h>
#include <float.h>

#include <string.h>
#include <stdlib.h>
// static long MAXINT = INT_MAX;

#define _GL_ON	1
#define _GL_OFF	0
#define _GL_NIL	0
#define PUBLIC
#define ABS(a)		(((a)>0)? (a):-(a))
#define MAXSTR	1024
#define NSD		3	/* Number of spatial dimensions */
#define MD(i,j)		((i)*NSD+(j))		/* 3-dim. Array with NSD */
#define M_LN(e,n)	((e)*8+(n))  /* for Lnods */
#define M_LN_n9(e,n)	((e)*9+(n))  /* for Lnods */
#define M_LN_q1q1(e,n)	((e)*16+(n))  /* for Lnods */
#define M_LN_sh4(e,n)	((e)*4+(n))  /* for Lnods */
#define _GL_COUT	'#'
#define _GL_ESC		'\\'
#define PI		3.14159265	



#define  MSHE	"_F.msh"

static int getNextLine(char *, FILE *);
static int makeLine(char *, FILE *);
static int countWord(char *);

static char  wsp[] = { ' ', '\t', _GL_NIL };

int
readANSYS_CDB(int &nNode, int &nElem, int* &Lnods, double* &InitGrid, std::string cdbFileName)
{
//	double *InitGrid;

	static FILE  *grid;
	int nElem_f, nElem_so;
	int /**Lnods,*/ *Lnods_so9;
	int *elem_n9;
	int  *Lnods_fq1q1, *Lnods_sh4, *elem_connect;
	double *grid_qun;
//	int nSupsDof_vx, nSupsDof_vy, nSupsDof_vz;
/*		nSupsDof_ux, nSupsDof_uy, nSupsDof_uz*/
/*		nSupsMesh_ux, nSupsMesh_uy, nSupsMesh_uz;*/
	char lab1[255], char2[255];
	int num1, nEelmTYPE, nMat;
	int *Node_finish;
	int i,/* nNode,*/ id,/* nElem,*/ ielem;
	int nodes8 = 8, nElem_sh;
	int idummy;
	double *Value_NodalForce;
	int izero, ione = 1;   
	char   str[MAXSTR];
	FILE  *fptemp, *temp1/* *fptemp_ref,*/ ;
	int  /*nSupsNode_fix_so,*/ *NodeType, *NodeType_qun, *ElemType,
		*Node_NodalForce, *DofIndex_NodalForce;
	int FluidEtype;
	int Q1Q1, Q1Q0, Q1;
	int MeshType = 7;

	Q1Q1 = 1;
	Q1Q0 = 2;
	Q1 = 3;



	FluidEtype = Q1;

	izero = 0;


	std::string workingDirection = "D:\OpenCASCADE-7.2.0-vc10-64\opencascade-7.2.0\samples\mfc\standard\win64";
	temp1 = fopen("temp", "w");//打开文件
	//grid = fopen(&(workingDirection + cdbFileName)[0]
	if ((grid = fopen("file.cdb", "r")) == NULL) {
//		printf("file open is failed");
//		return 0;
	}

	do {
		getNextLine(str, grid);
		sscanf(str, "%6s,%4s%*s%d", &lab1, &char2, &num1);
		if (!strcmp(lab1, "FINISH"))  exit(0);
		fprintf(temp1, " %s %s %d \n", lab1, char2, num1);
		if (!strcmp(char2, "NODE"))
		{
			nNode = num1;
		}
		if (!strcmp(char2, "ELEM"))
		{
			nElem = num1;
		}
		if (!strcmp(char2, "TYPE"))
		{
			nEelmTYPE = num1;
		}
		if (!strcmp(char2, "MAT"))
		{
			nMat = num1;
		}
	} while (strcmp(lab1, "NBLOCK") != 0);//NBLOCK之后为 具体的节点坐标信息，之前为控制信息

	printf("nNode %d\n", nNode);
	printf("nElem %d\n", nElem);

	NodeType = (int *)malloc((unsigned)(nNode)* sizeof(int));
	//	NodeType = (int *)malloc((unsigned)(nNode), sizeof(int));

	InitGrid = (double *)malloc((unsigned)(NSD*nNode) * sizeof(double));
	NodeType_qun = (int *)malloc((unsigned)(nNode * 3) * sizeof(int));
	grid_qun = (double *)malloc((unsigned)(NSD*nNode * 3) * sizeof(double));
	Node_finish = (int *)malloc((unsigned)(nNode)* sizeof(int));
	Node_NodalForce = (int *)malloc((unsigned)(nNode * 3) * sizeof(int));
	DofIndex_NodalForce = (int *)malloc((unsigned)(nNode * 3) * sizeof(int));
	Value_NodalForce = (double *)malloc((unsigned)(nNode * 3) * sizeof(double));

	Lnods = (int *)malloc((unsigned)(nElem * 8) * sizeof(int));
	ElemType = (int *)malloc((unsigned)(nElem)* sizeof(int));
	elem_connect = (int *)malloc((unsigned)(nElem)* sizeof(int));
	Lnods_fq1q1 = (int *)malloc((unsigned)(nElem * 16) * sizeof(int));
	Lnods_sh4 = (int *)malloc((unsigned)(nElem * 4) * sizeof(int));
	Lnods_so9 = (int *)malloc((unsigned)(nElem * 9) * sizeof(int));
	elem_n9 = (int *)malloc((unsigned)(nElem)* sizeof(int));

	getNextLine(str, grid);

	//由于cdb默认 0 坐标 输出为空，使得读入的这些点变为 个大数，需要处理，即，初始化initgrid为0列向量
	for (int i = 0; i < nNode; i++){
		for (int j = 0; j < 3;j++) InitGrid[MD(i, j)] = 0.0;		
	}


	for (i = 0; i<nNode; i++) {
		getNextLine(str, grid);
		sscanf(str, "%d %d %d %lf %lf %lf",
			&id, &idummy, &idummy,
			&InitGrid[MD(i, 0)], &InitGrid[MD(i, 1)], &InitGrid[MD(i, 2)]);
		fprintf(temp1, "%d %lf %lf %lf \n", i,
			InitGrid[MD(i, 0)], InitGrid[MD(i, 1)], InitGrid[MD(i, 2)]);
	}

	printf(" nElem %d \n", nElem);
	ielem = 0;
	nElem_f = 0;
	nElem_so = 0;
	nElem_sh = 0;

	do {
		getNextLine(str, grid);
		sscanf(str, "%6s,%s%*s%d", &lab1, &char2, &num1);
		/*    printf(" %s %s %d \n", lab1,char2, num1);*/
	} while (strcmp(lab1, "EBLOCK") != 0);

	getNextLine(str, grid);

	/*  Read connectivity */
	for (ielem = 0; ielem<nElem; ielem++) {
		getNextLine(str, grid);
		sscanf(str,
			"%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
			&idummy, &ElemType[ielem], &elem_connect[ielem], &idummy, &idummy,
			&idummy, &idummy, &idummy, &idummy, &idummy,
			&idummy,
			&Lnods[M_LN(ielem, 0)], &Lnods[M_LN(ielem, 1)],
			&Lnods[M_LN(ielem, 2)], &Lnods[M_LN(ielem, 3)],
			&Lnods[M_LN(ielem, 4)], &Lnods[M_LN(ielem, 5)],
			&Lnods[M_LN(ielem, 6)], &Lnods[M_LN(ielem, 7)]);

		if (ElemType[ielem] == 1)
		{
			nElem_f += 1;
		}
		else
		{
		}
	}

	// read components
	do {
		char seps[] = "\n";
		char seps1[] = ", \t";

		getNextLine(str, grid);
		sscanf(str, "%6s,%s%*s%d", &lab1, &char2, &num1);


	} while (strcmp(lab1, "FINISH") != 0);


	printf(" nElem_f %d ,nElem_so %d ,nElem_sh %d \n",
		nElem_f, nElem_so, nElem_sh);

	/* Open out put files */
	fptemp = fopen("meshin.msh", "w");


	/* Output grid */

	fprintf(fptemp, "%d   %d   %d   %d \n", NSD, nNode, nElem, nodes8);
	for (i = 0; i<nNode; i++) {
		fprintf(fptemp, "           %d    %.5e    %.5e    %.5e \n", i + 1,
			InitGrid[MD(i, 0)], InitGrid[MD(i, 1)],
			InitGrid[MD(i, 2)]);
	}


	/* Output element connectivities */



	for (ielem = 0; ielem<nElem; ielem++) {

		fprintf(fptemp, "            %d   %d   %d         %d    %d    %d    %d    %d    %d    %d    %d \n",
			ielem + 1, nodes8, MeshType,
			Lnods[M_LN(ielem, 0)], Lnods[M_LN(ielem, 1)],
			Lnods[M_LN(ielem, 2)], Lnods[M_LN(ielem, 3)],
			Lnods[M_LN(ielem, 4)], Lnods[M_LN(ielem, 5)],
			Lnods[M_LN(ielem, 6)], Lnods[M_LN(ielem, 7)]);
	}

	fprintf(fptemp, " %d \n ", ione);

	fprintf(fptemp, " Displacement  \n ");
	for (i = 0; i<nNode; i++) {
		fprintf(fptemp, "  %.5e    \n", (double)i);
	}

	fclose(fptemp);

//	system("pause");

	return 1;
}

/*
*  getline.c  --- Read uncommented one line, and set strings
*
*	Written by kawakami@sml.t.u-tokyo.ac.jp
*	$Id: getline.c,v 0.3 1998-10-06 11:38:58+09 kawakami Exp $
*/

static int
getNextLine(char  *str, FILE  *fp)			/* P_( (char *,FILE *) ) */
{
	int  ret;

	if (makeLine(str, fp))
		ret = countWord(str);
	else
		ret = -1;
	return(ret);
}


static int
makeLine(char  *str, FILE  *fp)			/* P_( (char *,FILE *) ) */
{
	char  line[MAXSTR], *strp;
	int   i;

	strp = str;
	for (;;) {
		if (fgets(line, MAXSTR, fp) == NULL)
			return(0);
		if (line[0] == _GL_COUT || line[0] == '\n')
			continue;
		i = 0;
		while (i<MAXSTR) {
			if (line[i] == _GL_ESC) {
				i++;
				switch (line[i]) {
				case _GL_COUT:
					*strp++ = _GL_COUT;
					break;
				case _GL_NIL:
					*strp = _GL_NIL;
					return(1);
				case _GL_ESC:
					*strp++ = _GL_ESC;
					break;
				case '\n':
					goto next;
				default:
					*strp++ = _GL_ESC;
					*strp++ = line[i];
				}
			} /* if line[i] */
			else if (line[i] == '\n' || line[i] == _GL_NIL || line[i] == _GL_COUT) {
				*strp = _GL_NIL;
				return(1);
			}
			else
				*strp++ = line[i];
			i++;
		} /* while i */
	next:
		;
	}	/* for */
}


static  int
countWord(char  *str)			/* P_( (char *) ) */
{
	int  word, word_on, i, j;

	word = 0;
	word_on = _GL_OFF;
	i = 0;
	while (str[i] != _GL_NIL) {
		j = 0;
		while (wsp[j] != _GL_NIL) {
			if (str[i] == wsp[j]) {
				word_on = _GL_OFF;
				goto next;
			}
			j++;
		}
		if (!word_on) {
			word_on = _GL_ON;
			word++;
		}
	next:
		i++;
	} /* while str[i] */
	return(word);
}

/* end of file */
