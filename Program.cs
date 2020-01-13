using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;


namespace watgay
{
	public static class gay
	{
	
		public static bool cplxmat(string input, string pattern)
		{

			int[] inputPosStack = new int[(input.Length + 1) * (pattern.Length + 1)];   // Stack containing input positions that should be tested for further matching
			int[] patternPosStack = new int[inputPosStack.Length];                      // Stack containing pattern positions that should be tested for further matching
			int stackPos = -1;                                                          // Points to last occupied entry in stack; -1 indicates that stack is empty
			bool[,] pointTested = new bool[input.Length + 1, pattern.Length + 1];        // Each true value indicates that input position vs. pattern position has been tested

			int inputPos = 0;   // Position in input matched up to the first multiple wildcard in pattern
			int patternPos = 0; // Position in pattern matched up to the first multiple wildcard in pattern

			// Match beginning of the string until first multiple wildcard in pattern
			while (inputPos < input.Length && patternPos < pattern.Length && pattern[patternPos] != '*' && (input[inputPos] == pattern[patternPos] || pattern[patternPos] == '?')) {
				inputPos++;
				patternPos++;
			}

			// Push this position to stack if it points to end of pattern or to a general wildcard
			if (patternPos == pattern.Length || pattern[patternPos] == '*') {
				pointTested[inputPos, patternPos] = true;
				inputPosStack[++stackPos] = inputPos;
				patternPosStack[stackPos] = patternPos;
			}
			bool matched = false;

			// Repeat matching until either string is matched against the pattern or no more parts remain on stack to test
			while (stackPos >= 0 && !matched) {

				inputPos = inputPosStack[stackPos];         // Pop input and pattern positions from stack
				patternPos = patternPosStack[stackPos--];   // Matching will succeed if rest of the input string matches rest of the pattern

				if (inputPos == input.Length && patternPos == pattern.Length)
					matched = true;     // Reached end of both pattern and input string, hence matching is successful
        else {   
					// First character in next pattern block is guaranteed to be multiple wildcard
					// So skip it and search for all matches in value string until next multiple wildcard character is reached in pattern

					for (int curInputStart = inputPos; curInputStart < input.Length; curInputStart++) {

						int curInputPos = curInputStart;
						int curPatternPos = patternPos + 1;

						if (curPatternPos == pattern.Length) {   // Pattern ends with multiple wildcard, hence rest of the input string is matched with that character
							curInputPos = input.Length;
						} else {

							while (curInputPos < input.Length && curPatternPos < pattern.Length && pattern[curPatternPos] != '*' &&
							       (input[curInputPos] == pattern[curPatternPos] || pattern[curPatternPos] == '?')) {
								curInputPos++;
								curPatternPos++;
							}

						}

						// If we have reached next multiple wildcard character in pattern without breaking the matching sequence, then we have another candidate for full match
						// This candidate should be pushed to stack for further processing
						// At the same time, pair (input position, pattern position) will be marked as tested, so that it will not be pushed to stack later again
						if (((curPatternPos == pattern.Length && curInputPos == input.Length) || (curPatternPos < pattern.Length && pattern[curPatternPos] == '*'))
						    && !pointTested[curInputPos, curPatternPos]) {
							pointTested[curInputPos, curPatternPos] = true;
							inputPosStack[++stackPos] = curInputPos;
							patternPosStack[stackPos] = curPatternPos;
						}

					}
				}

			}

			return matched;

		}


		public static HTSstm[] prvtreez = new HTSstm[256];
		public static int nprvtreez = 0;

		public static int matzn = 0;
		public static matz[] mmt;
		public static int matzlim = 1102;

		public static class vindov
		{
			public static double[,] coeff = new double[3, 3] {
				{ 999, 1.0, 999 },
				{ -0.5, 0.0, 0.5 },
				{ 1.0, -2.0, 1.0 }
			};
			public static int[] winL = new int[3] {
				0,
				-1,
				-1
			};
			public static int[] winR = new int[3] {
				0,
				1,
				1
			};
			
		}
		

                                        
		public class sstream
		{
                                                            
			public double pitch = 0;
			public int dur;
			
			public ushort[] durseg;
			public double[][][,,] pamz = new double[10][][,,];
			public double[][] msdz = new double[10][];
			
			
			public void adjDur(int dest, double rho)
			{
				double temp1 = 0;
				double temp2 = 0;
				int j = 0;
				int sttlen = pamz[0][0].Length / 2;
				int lablen = pamz[0].Length;
				while (dur != dest) {
					if (dest > dur) {
						j = -1;
						for (int i = 0; i < lablen; i++) {
							var homo = pamz[0][i];
							for (int i2 = 0; i2 < sttlen; i2++) {
								temp2 = Math.Abs(rho - ((double)durseg[i] + 1 - homo[0, 0, i2]) / homo[1, 0, i2]);
								if (j < 0 || temp1 > temp2) {
									j = i * sttlen + i2;
									temp1 = temp2;
								}
							}
						}
						dur++;
						durseg[j]++;
					} else {
						j = -1;
						for (int i = 0; i < lablen; i++) {
							var homo = pamz[0][i];
							for (int i2 = 0; i2 < sttlen; i2++) {
								if (durseg[i] > 1) {
									temp2 = Math.Abs(rho - ((double)durseg[i] - 1 - homo[0, 0, i2]) / homo[1, 0, i2]);
									if (j < 0 || temp1 > temp2) {
										j = i * sttlen + i2;
										temp1 = temp2;
									}
								}
							}
						}
						dur--;
						durseg[j]--;
					}
				}
				
				
				return;
			}
                                                            
			public void setspeed(double speed)
			{
				int tstate = pamz[3].Length;
				durseg = new ushort[tstate];
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
				int sttlen = pamz[0][0].Length / 2;
				int lablen = pamz[0].Length;
				double durf = 0;                                                                
                                                                                
                                                                                
				if (speed == 1.0) {
					dur = 0;
					for (int i = 0; i < lablen; i++) {
						var homo = pamz[0][i];
						for (int j = 0; j < sttlen; j++) {
							durf += homo[0, 0, j];
							double t1 = homo[0, 0, j] + 0.5;
							short t1n = (short)t1;
							if (t1n < 1)
								t1n = 1;
                                                                                                                        
                            
							//=======================
							//if ((i * sttlen + j) == 31)
							//	t1n++;
							
							//=======================
							durseg[i * sttlen + j] = (ushort)t1n;
							dur += t1n;
                                                                                                                                 
                                                                                                                        
                                                                                                                        
						}
					}
					
					int idurf = (int)(durf + 0.5);
					
					if (dur != idurf) {
						adjDur(idurf, 0);
						return;
					}
					
					return;
				}
                                                                                
                                                                                
				
				double origvari = 0;
				foreach (var homo in pamz[0]) {
					for (int i = 0; i < sttlen; i++) {
						durf += homo[0, 0, i];
						origvari += homo[1, 0, i];
					}
				}
                                                                                
				double origmean = durf;
				durf /= speed;
                                                                                
				dur = (int)(durf + 0.5);
				if (dur < 1)
					dur = 1;
                                                                                                    
				if (dur < tstate) {
					for (int i = 0; i < tstate; i++)
						durseg[i] = 1;
                                                                                                                        
                                                                                                                        
					dur = tstate;
					return;
				}
                                                                               
				durf = (durf - origmean) / origvari; //rho=tmpdurf
                                                                               
				int sum = 0;
				for (int i = 0; i < lablen; i++) {
					var homo = pamz[0][i];
					for (int j = 0; j < sttlen; j++) {
						double t1 = homo[0, 0, j] + durf * homo[1, 0, j] + 0.5;
						short t1n = (short)t1;
						if (t1n < 1)
							t1n = 1;
                                                                                                                        
						durseg[i * sttlen + j] = (ushort)t1n;
						sum += t1n;
                                                                                                                        
                                                                                                                        
					}
				}
				
				
				if (sum != dur) {
					int idur = dur;
					dur = sum;
					adjDur(idur, durf);
					return;
				}
				
                                                                                
				/*
                                                                                while(sum!=dur)
                                                                                {
                                                                                                   Console.WriteLine(dur+" but "+sum);
                                                                                }
                                                                                */
			}
			
			public sstream(HTSvoice[] voiz, labio[] labz)
			{
				int shuaa = labz.Length;
				int nvoi = voiz.Length;
			
			
					
				
				for (int j = 0; j < 10; j++) {

					if (voiz[0].stmz[j] == null)
						continue;

                                                                                                    
					int ntrr = voiz[0].stmz[j].HTStr.Length;
					pamz[j] = new double[shuaa * ntrr][,,];
                                                                                                   
					bool ismsd = voiz[0].stmz[j].msd;
					if (ismsd)
						msdz[j] = new double[shuaa * ntrr];
					for (int vvo = 0; vvo < nvoi; vvo++) {
						//for(int i=0;i<shuaa;i++)
						Parallel.For(0, shuaa, i => {
					
							for (int k = 0; k < ntrr; k++) {
							//Parallel.For(0, ntrr, k => {
								if (ismsd)
									pamz[j][ntrr * i + k] = voiz[vvo].stmz[j].domati(labz[i], k, pamz[j][ntrr * i + k], ref msdz[j][ntrr * i + k], 1.0);
								else {
									double shybm = 0;
									pamz[j][ntrr * i + k] = voiz[vvo].stmz[j].domati(labz[i], k, pamz[j][ntrr * i + k], ref shybm, 1.0);
								}
                                                                                                                                            
						
							}//);
						}); 
					}
                                                            
               
					if (j > 2) {
						int allfrm = ntrr * shuaa;
						Parallel.For(0, allfrm, i => {
				             	
							double[,,] myke = pamz[j][i];
							int vecll = myke.GetLength(2);
							//Parallel.For(0, 3, win => {
							for (int win = 0; win < 3; win++) {
								Parallel.For(0, vecll, vec => {
									double x = myke[1, win, vec];
									double ot = 0;
				             	             	             	
									if (x >= 1.0e+19)
										ot = 0.0;
									else if (x <= -1.0e+19)
										ot = 0.0;
									else if (x <= 1.0e-19 && x >= 0)
										ot = 1.0e+38;
									else if (x >= -1.0e-19 && x < 0)
										ot = -1.0e+38;
									else
										ot = (1.0 / x);
				             	             	             	
									myke[1, win, vec] = ot;
				             	             	             	
								});
				             	             	
				             	             	
							}//);
				             	
				             	
						});
					}
                                                            
				}
			}
			
			
			double[,,] tysfo(double[,,] inin, int veclr)
			{
				double[,,] nyuss = (double[,,])inin.Clone();
				for (int k = 1; k < 3; k++) {
									
					for (int ll = 0; ll < veclr; ll++)
						nyuss[1, k, ll] = 0;
									
											
				}
			
				return nyuss;
			}
			
			public double[][,] produce()
			{
				int[] ydax = new int[5];
				int npstm = 0;
			
				for (int i = 3; i < 10; i++) {
					if (pamz[i] != null) {
						ydax[npstm] = i;
						npstm++;
					}
                                                                                                    
				}
			 
				double[][,] gstrm = new double[npstm][,];
			
			
				//for (int i = 0; i < npstm; i++) {
				Parallel.For(0, npstm, i => {
			             	
					double[] msd = msdz[ydax[i]];
					double[][,,] nauss = pamz[ydax[i]];
					int veclr = nauss[0].GetLength(2);
				
					int tstate = nauss.Length;
					int tmpdur = dur;
					gstrm[i] = new double[dur, veclr];
					bool[] tmsdfz = null;
				
				
				
					if (msd != null) {
						tmsdfz = new bool[dur];
					
						tmpdur = 0;
						int naofrm = 0;
						for (int i1 = 0; i1 < tstate; i1++) {
						
							int naodurseg = durseg[i1];
						
						
							if (msd[i1] > 0.5) {
								tmpdur += naodurseg;
							
								for (int i2 = 0; i2 < naodurseg; i2++) {
									tmsdfz[naofrm + i2] = true;
								}
							
							} else {
								for (int i2 = 0; i2 < naodurseg; i2++) {
									tmsdfz[naofrm + i2] = false;
								}
							}
						
							naofrm += naodurseg;
							
						}
					
					
					}
				
					double[][,,] ptmp = new double[tmpdur][,,];
				
					if (msd != null) {
						int shortm = 0;
						int naofrm = 0;
						for (int i1 = 0; i1 < tstate; i1++) {
							int naodurseg = durseg[i1];
							if (tmsdfz[naofrm]) {
							
								int sta = 0;
								if (naofrm == 0 || !tmsdfz[naofrm - 1]) {
								
									ptmp[shortm] = tysfo(nauss[i1], veclr); //tysfo(,veclr);
									sta = 1;
								}
							
								for (int i2 = sta; i2 < naodurseg; i2++) {
									/*
								bool some0 = false;
								bool[] set0 = new bool[3];
								
								for (int k = 0; k < 3; k++) {
									bool not_bound = true;
									for (int shift = gay.vindov.winL[k]; shift <= gay.vindov.winR[k]; shift++) {
										if ((naofrm + i2 + shift) < 0 || sstmj.dur <= (naofrm + i2 + shift) || !tmsdfz[naofrm + i2 + shift]) {
											not_bound = false;
											break;
										}
									}
									
									if (!not_bound && k != 0) {
										some0 = true;
										set0[k] = true;
									}
									
								}
								
								if (some0) {
									double[,,] nyuss = (double[,,])nauss[i1].Clone();
									for (int k = 0; k < 3; k++) {
										if (set0[k]) {
											for (int ll = 0; ll < veclr; ll++)
												nyuss[1, k, ll] = 0;
										}
											
									}
									ptmp[shortm + i2] = nyuss;
								} else {
									*/
									ptmp[shortm + i2] = nauss[i1];
									//}
								
								}
							
							
							
								shortm += naodurseg;
							} else {
								
								ptmp[shortm - 1] = tysfo(ptmp[shortm - 1], veclr);
							
							}
					
						
							naofrm += naodurseg;
						}
					
					
						ptmp[tmpdur - 1] = tysfo(ptmp[tmpdur - 1], veclr);
					
					} else {
						int naofrm = 0;
						for (int i1 = 0; i1 < tstate; i1++) {
							int naodurseg = durseg[i1];
							for (int i2 = 0; i2 < naodurseg; i2++) {
								/*
							bool some0 = false;
							bool[] set0 = new bool[3];
							
							
							for (int k = 0; k < 3; k++) {
									bool not_bound = true;
									for (int shift = gay.vindov.winL[k]; shift <= gay.vindov.winR[k]; shift++) {
										if ((naofrm + i2 + shift) < 0 || sstmj.dur <= (naofrm + i2 + shift)) {
											not_bound = false;
											break;
										}
									}
									
									if (!not_bound && k != 0) {
										some0 = true;
										set0[k] = true;
									}
									
								}
							
							
							if (some0) {
								double[,,] nyuss = (double[,,])nauss[i1].Clone();
								for (int k = 0; k < 3; k++) {
									if (set0[k]) {
										for (int ll = 0; ll < veclr; ll++)
											nyuss[1, k, ll] = 0;
									}
											
								}
								ptmp[naofrm + i2] = nyuss;
							} else {
								*/
								ptmp[naofrm + i2] = nauss[i1];
								//}
							
							}
					
							naofrm += naodurseg;
						}
					
					
						for (int i1 = 0; i1 < tmpdur; i1 += (tmpdur - 1)) {
								
							ptmp[i1] = tysfo(ptmp[i1], veclr);
					
						}
					
					}
                  
					//=====================
				
					/*
				for (int i1 = 0; i1 < tmpdur; i1++) {
					for (int i2 = 0; i2 < 3; i2++) {
						for (int i3 = 0; i3 < veclr; i3++) {
							Console.WriteLine("stm" + i + ", frm" + i1 + ", win" + i2 + ", vec" + i3 + ", mean= " + ptmp[i1][0, i2, i3].ToString("#0.000000"));
							Console.WriteLine("stm" + i + ", frm" + i1 + ", win" + i2 + ", vec" + i3 + ", vari= " + ptmp[i1][1, i2, i3].ToString("#0.000000"));
						}
					}
				}
				*/
				

					//=====================
				
					double[,] paris = new double[tmpdur, veclr];
				
					Parallel.For(0, veclr, vecN => {
						//for (int vecN = 0; vecN < veclr; vecN++) {
					
						double[] wum = new double[tmpdur];
						double[,] wuw = new double[tmpdur, 3];
					
						for (int t = 0; t < tmpdur; t++) {
							for (int winN = 0; winN < 3; winN++) {
								for (int shift = vindov.winL[winN]; shift <= vindov.winR[winN]; shift++) {
									if ((t + shift >= 0) && (t + shift < tmpdur) && (vindov.coeff[winN, 1 - shift] != 0.0)) {
									
										double wu = vindov.coeff[winN, 1 - shift] * ptmp[t + shift][1, winN, vecN];
										wum[t] += wu * ptmp[t + shift][0, winN, vecN];
               						
										for (int jj = 0; (jj < 3) && (t + jj < tmpdur); jj++)
											if ((jj <= vindov.winR[winN] + shift) && (vindov.coeff[winN, 1 + jj - shift] != 0.0))
												wuw[t, jj] += wu * vindov.coeff[winN, 1 + jj - shift];
									}
								}
							}
						}
					
					
						for (int t = 0; t < tmpdur; t++) {
							for (int winN = 1; (winN < 3) && (t >= winN); winN++) {
								wuw[t, 0] -= wuw[t - winN, winN] * wuw[t - winN, winN] * wuw[t - winN, 0];
							
							}
						
							for (int winN = 1; winN < 3; winN++) {
								for (int jj = 1; (winN + jj < 3) && (t >= jj); jj++)
									wuw[t, winN] -= wuw[t - jj, jj] * wuw[t - jj, winN + jj] * wuw[t - jj, 0];
							
								wuw[t, winN] /= wuw[t, 0];
							}
						
						}
					
						double[] g = new double[tmpdur];
					
						for (int t = 0; t < tmpdur; t++) {
							g[t] = wum[t];
							for (int winN = 1; (winN < 3) && (t >= winN); winN++)
								g[t] -= wuw[t - winN, winN] * g[t - winN];
						}
					
						for (int rev = 0; rev < tmpdur; rev++) {
							int t = tmpdur - 1 - rev;
							paris[t, vecN] = g[t] / wuw[t, 0];
							for (int winN = 1; (winN < 3) && (t + winN < tmpdur); winN++)
								paris[t, vecN] -= wuw[t, winN] * paris[t + winN, vecN];
						}
					
					});
				
					//=======================
				
					/*
					for (int t = 0; t < tmpdur; t++) {
						string siya = "stm" + i + ", t" + t + ": " + '\t';
					
						for (int vv = 0; vv < veclr; vv++) {
							siya += paris[t, vv].ToString("#0.000000") + ", " + '\t';
						}
					
						Console.WriteLine(siya);
					
					}
					*/
				
					//=======================
					
					if (msd != null) {
						double[,] rpar = new double[dur, veclr];
						int msdfrm = 0;
						for (int i1 = 0; i1 < dur; i1++) {
							if (tmsdfz[i1]) {
								for (int i2 = 0; i2 < veclr; i2++) {
									rpar[i1, i2] = paris[msdfrm, i2];
								}
								msdfrm++;
							} else {
								for (int i2 = 0; i2 < veclr; i2++) {
									rpar[i1, i2] = -1.0e+10;
								}
							}
						}
						
						gstrm[i] = rpar;
						paris = new double[0, 0];
					} else
						gstrm[i] = paris;
				
				});
				
				
				return gstrm;
			}
                                                            
		}
                                        
                                        
  
		public class labio
		{
			public string lab;
			public byte[] matpi;
                                                            
			public labio(string labr)
			{
				lab = labr;
				matpi = new byte[gay.matzlim];
				matpi[0] = 0xff;
				matpi[1] = 1;
			}
		}


		public class HTSvoice
		{
			public int statesz;
			public HTSstm[] stmz;
                    
                    
			internal void readPDUR(string bdir, string inty, int stmn)
			{
				stmz[stmn] = new HTSstm(inty, 1, 1, 2 * 4, bdir, false);
				int[] ofz = new int[1];
				ofz[0] = 4;
				stmz[stmn].qsus(ofz);

			}
                    
			public HTSvoice(string bdir)
			{
				bdir += '\\';
				statesz = 5;
				stmz = new HTSstm[10];
                                        
				stmz[0] = new HTSstm("DUR", statesz, 1, 2 * statesz * 4, bdir, false);
				int[] ofz = new int[1];
				ofz[0] = 4;
				stmz[0].qsus(ofz);
                                        
				if (File.Exists(bdir + "pdr.pd"))
					readPDUR(bdir, "PDR", 1);

				if (File.Exists(bdir + "rlc.pd"))
					readPDUR(bdir, "RLc", 1);

				if (File.Exists(bdir + "rls.pd"))
					readPDUR(bdir, "RLs", 2);
					
					
				bool msd = false;
				int vecn = 50;
				int winn = 3;
				int pdsz = (vecn * winn * 2) * 4;

				stmz[3] = new HTSstm("MGC", vecn, winn, pdsz, bdir, msd);


				ofz = new int[statesz];
				ofz[0] = statesz * 4;

				for (int i = 0; i < statesz - 1; i++) {
					ofz[i + 1] = ofz[i] + stmz[3].noka(i * 4) * pdsz;
				}

				stmz[3].qsus(ofz);
				//============
                                                                                
				msd = true;
				winn = 3;
				vecn = 1;
				pdsz = (vecn * winn * 2 + 1) * 4;
				stmz[4] = new HTSstm("LF0", vecn, winn, pdsz, bdir, msd);


				ofz = new int[statesz];
				ofz[0] = statesz * 4;

				for (int i = 0; i < statesz - 1; i++) {
					ofz[i + 1] = ofz[i] + stmz[4].noka(i * 4) * pdsz;
				}
				stmz[4].qsus(ofz);
				//============
   
				msd = false;
				winn = 3;
				vecn = 22;
				pdsz = vecn * winn * 2 * 4;
				stmz[5] = new HTSstm("BAP", vecn, winn, pdsz, bdir, msd);


				ofz = new int[statesz];
				ofz[0] = statesz * 4;

				for (int i = 0; i < statesz - 1; i++) {
					ofz[i + 1] = ofz[i] + stmz[5].noka(i * 4) * pdsz;
				}

				stmz[5].qsus(ofz);
				//============
   
				if (File.Exists(bdir + "vib.pd")) {
					msd = true;
					winn = 3;
					vecn = 2;
					pdsz = (vecn * winn * 2 + 1) * 4;

					stmz[6] = new HTSstm("VIB", vecn, winn, pdsz, bdir, msd);


					ofz = new int[statesz];
					ofz[0] = statesz * 4;

					for (int i = 0; i < statesz - 1; i++) {
						ofz[i + 1] = ofz[i] + stmz[6].noka(i * 4) * pdsz;
					}

					stmz[6].qsus(ofz);
   
				}
                                        
                                        
			}
		}


		public class matz
		{
			public string mya;
			public bool isSimp;
			public bool fromhead;
			public bool fromtail;
			public bool hasq;
			public bool hasX;

		}

		public class HTSstm
		{

			public string stmsig;
			public brnchy[][] HTStr;
			public byte[] pdbuf;
			public byte[] aly;
			public bool msd = false;
   
			public int vecn;
			public int winn;
			public int pdsz;
			public keyta[] lyto;
			public matz[] mmts;
			public ulong hashy;


			public HTSstm(string syg, int vecnr, int winnr, int pdszr, string bdir, bool msdr)
			{
				msd = msdr;
				stmsig = syg;
				vecn = vecnr;
				winn = winnr;
				pdsz = pdszr;
				mmts = gay.mmt;

				pdbuf = File.ReadAllBytes(bdir + stmsig + ".pd");
				aly = File.ReadAllBytes(bdir + stmsig + ".tr.biz");
				mkhash();


			}

			unsafe void mkhash()
			{
				ulong hahaha = 0;
				int lyp = (aly.Length >> 3);
				fixed(byte * nrf = &(aly[0])) {
					for (int i = 0; i < lyp; i++) {
						ulong auu = *((ulong*)&(nrf[i * 8]));
						hahaha ^= auu;
					}


				}
				hashy = hahaha;

			}

			public unsafe int noka(int pyz)
			{
				fixed(byte * nrf = &(pdbuf[pyz])) {
					return *((int*)&(nrf[0]));
				}
			}

			public brnchy creifnul(brnchy inin)
			{
				if (inin == null)
					return new brnchy();
				else
					return inin;
			}

			public unsafe int treeing(int curpy, int ntrr, int skpfirst)
			{

				ushort ntreez;
				fixed(byte * nrf = & (aly[curpy])) {
					ntreez = *((ushort*)&(nrf[0]));
				}
				curpy += 2;

				HTStr[ntrr] = new brnchy[ntreez];
				brnchy[] otri = HTStr[ntrr];



				fixed(byte * nrf = & (aly[curpy])) {
					for (int i = 0; i < ntreez; i++) {
						ushort whq = *((ushort*)&(nrf[i * 6]));
						short noq = *((short*)&(nrf[i * 6 + 2]));
						short yesq = *((short*)&(nrf[i * 6 + 4]));

						otri[i] = creifnul(otri[i]);
						otri[i].key = lyto[whq];
						if (noq < 0) {
							int tadi = -noq;
							otri[tadi] = creifnul(otri[tadi]);
							otri[i].ifn = otri[tadi];
						} else {
							otri[i].pdn = skpfirst + noq * pdsz;
						}

						if (yesq < 0) {
							int tadi = -yesq;
							otri[tadi] = creifnul(otri[tadi]);
							otri[i].ify = otri[tadi];
						} else {
							otri[i].pdy = skpfirst + yesq * pdsz;
						}

					}
				}

				curpy += (6 * ntreez);

				//sinope(otri[0], "T" + ntrr);
   
				return curpy;



			}



			public unsafe void qsus(int[] ofz)
			{
  
  
  
				for (int i = 0; i < nprvtreez; i++) {
					if (hashy == prvtreez[i].hashy) {
						lyto = prvtreez[i].lyto;
						HTStr = prvtreez[i].HTStr;
						aly = new byte[0];
						return;
					}
				}
  
				prvtreez[nprvtreez] = this;
				nprvtreez++;
  

				ushort nq;
				fixed(byte * nrf = & (aly[0])) {
					nq = *((ushort*)&(nrf[0]));
				}

				lyto = new keyta[nq];

				int curpy = 0;
				fixed(byte * nrf = & (aly[2])) {

					for (int i = 0; i < nq; i++) {
						byte nmatz = *((byte*)&(nrf[curpy]));
						keyta hayi = new keyta();
     
     
						//hayi.tit = "Ques" + i; //cmtoutthis
     
     
						hayi.mati = new ushort[nmatz];
						curpy += 1;

						for (int ppy = 0; ppy < nmatz; ppy++) {
							hayi.mati[ppy] = *((ushort*)&(nrf[curpy + 2 * ppy]));
						}
						lyto[i] = hayi;

						curpy += (2 * nmatz);

					}
				}

				//Console.WriteLine("nupo: "+(curpy+2).ToString("X8"));
				//Console.ReadKey();

				curpy += 2;


				int lyppa = ofz.Length;
				HTStr = new brnchy[lyppa][];


				for (int i = 0; i < lyppa; i++) {

					curpy = treeing(curpy, i, ofz[i]);


				}
   
				aly = new byte[0];


			}

			void loggi(string patan, string typi, string labi)
			{
				Console.WriteLine(patan + "是" + labi);
			}
                                                            
			bool yRn(labio inlab, brnchy inbrn)
			{
				bool noNext = true;
				foreach (ushort mk in  inbrn.key.mati) {
					if (inlab.matpi[mk] == 1)
						return true;
					else if (inlab.matpi[mk] == 0)
						noNext = false;
				}
				
				if (noNext)
					return false;

				foreach (ushort mk in  inbrn.key.mati) {
				
					if (inlab.matpi[mk] == 1)
						return true;
				
					if (inlab.matpi[mk] == 0) {
						matz mtnow = mmts[mk];
						if (mtnow.isSimp) {
							if (inlab.lab.Contains(mtnow.mya)) {
								inlab.matpi[mk] = 1;
								//loggi('*'+mtnow.mya+'*',"anyplace",inlab.lab);
								return true;
							} else {
								inlab.matpi[mk] = 0xFF;
							}
						} else if (!mtnow.hasq && !mtnow.hasX) { 
							if (mtnow.fromhead) {
								if (inlab.lab.StartsWith(mtnow.mya)) {
									inlab.matpi[mk] = 1;
									//loggi(mtnow.mya+'*',"start",inlab.lab);
									return true;
								} else {
									inlab.matpi[mk] = 0xFF;
								}
							} else if (mtnow.fromtail) {
								if (inlab.lab.EndsWith(mtnow.mya)) {
									inlab.matpi[mk] = 1;
									//loggi('*'+mtnow.mya,"xx",inlab.lab);
									return true;
								} else {
									inlab.matpi[mk] = 0xFF;
								}
							} else {
								Console.WriteLine("wtf?" + mtnow.mya);
								Console.ReadKey();
							}

						} else {
							//complexhere
							string patan = mtnow.mya;
							if (mtnow.fromhead)
								patan += '*';
							else if (mtnow.fromtail)
								patan = '*' + patan;
							else
								patan = '*' + patan + '*';
							
							if (gay.cplxmat(inlab.lab, patan)) {
								inlab.matpi[mk] = 1;
								//loggi(patan,"cmplx",inlab.lab);
								return true;
							} else {
								inlab.matpi[mk] = 0xFF;
							}
						}


					}
										
				}

				return false;

			}


			int sashi(labio inlab, brnchy inbrn)
			{
				if (yRn(inlab, inbrn)) {
					if (inbrn.pdy == -1)
						return sashi(inlab, inbrn.ify);
					else
						return inbrn.pdy;
				} else {
					if (inbrn.pdn == -1)
						return sashi(inlab, inbrn.ifn);
					else
						return inbrn.pdn;


				}
				
				
				


			}

			//bool needsib=true;
			public unsafe double[,,] domati(labio inlab, int Tn, double[,,] prev, ref double msdp, double wegt)
			{
				int pdx = sashi(inlab, HTStr[Tn][0]);
				
				//=========
				/*
				int myir=0;
				if(Tn>0)
				{
                                                                                                    fixed(byte* pio = &(pdbuf[0]))
                                                                                                    {
                                                                                                    
                                                                                                    int* shka = (int*)pio;
                                                                                                   if(needsib)
                                                                                                   { 
                                                                                                    for(int i=1;i<5;i++)
                                                                                                    {
                                                                                                                        shka[i]+=shka[i-1];
                                                                                                                        
                                                                                                                        //Console.WriteLine("shk: "+shka[i]);
                                                                                                    }
                                                                                                    needsib=false;
                                                                                                    }
                                                                                                    myir=shka[Tn-1];
                                                                                                    }
                                                                                                    
                                                                                                    
				}
				
				myir=(pdx/pdsz)-myir+1;
				
				
				Console.WriteLine(this.stmsig+'_'+Tn+'='+myir);
				*/
				//=========
				
				
				
				if (prev == null) {
					prev = new double[2, winn, vecn];
				}
				
				fixed(byte* pio = &(pdbuf[pdx])) {
					float* shka = (float*)pio;
                                                                                                    
                                                                                                    
                                                                                                    
					for (int aa = 0; aa < 2; aa++) {
						for (int bb = 0; bb < winn; bb++) {
							for (int cc = 0; cc < vecn; cc++) {
								prev[aa, bb, cc] += shka[(aa * winn + bb) * vecn + cc] * wegt;
							}    
						}
					}
                                                                                                                       
                                                                                                                        
                                                                                                                        
					if (msd)
						msdp = shka[2 * vecn * winn] * wegt;
				}
				
				return prev;
				
			}

		}


		public class brnchy
		{
			public keyta key;
			public brnchy ify;
			public brnchy ifn;
			public int pdy = -1;
			public int pdn = -1;
		}

		public class keyta
		{
			public string tit;
			public ushort[] mati;



			public void makefol(string pathh)
			{
				Directory.CreateDirectory(pathh);
				int nsigua = mati.Length;
				string[] sigua = new string[nsigua];
				for (int i = 0; i < nsigua; i++) {
					matz tmptz = gay.mmt[mati[i]];
					string ampi = tmptz.mya;
					if (tmptz.isSimp) {
						ampi = '*' + ampi + '*';
					} else {
						if (tmptz.fromhead) {
							ampi = ampi + '*';
						} else if (tmptz.fromtail) {
							ampi = '*' + ampi;
						} else {
							ampi = '*' + ampi + '*';
						}

					}



					sigua[i] = ampi;
				}


				File.WriteAllLines(pathh + '\\' + tit + ".q", sigua);

			}
		}



		static void rydprv()
		{

			mmt = new matz[matzlim];

			mmt[0] = new matz(); //alwaysfalse


			mmt[1] = new matz(); //alwaysture


			matzn = 2;

			byte[] saite = File.ReadAllBytes("mtzprv.sav");
			int lyn = matzlim - 2; //saite[0]+saite[1]*256;
			int nypo = 2;
			for (int i = 0; i < lyn; i++) {
				mmt[2 + i] = new matz();
				int urai = saite[nypo];
				byte sig = (byte)(urai & 0xf0);
				urai = urai & 0xF;

				if (sig == 0)
					mmt[2 + i].isSimp = true;
				else {
					mmt[2 + i].isSimp = false;
					if ((sig & 0x10) != 0) {
						if ((sig & 0x20) != 0)
							mmt[2 + i].fromhead = true;
						else
							mmt[2 + i].fromtail = true;
					} else {
						mmt[2 + i].fromhead = false;
						mmt[2 + i].fromtail = false;
					}

					if ((sig & 0x40) != 0) {
						if ((sig & 0x80) != 0) {
							mmt[2 + i].hasX = true;
							mmt[2 + i].hasq = false;
							for (int kk = 0; kk < urai; kk++) {
								if (saite[nypo + 1 + kk] == '?') {
									mmt[2 + i].hasq = true;
									break;
								}


							}
						} else {
							mmt[2 + i].hasq = true;
							mmt[2 + i].hasX = false;
						}
					}




				}


				mmt[2 + i].mya = System.Text.Encoding.ASCII.GetString(saite, nypo + 1, urai);
				nypo += (1 + urai);
			}

			matzn += lyn;


		}



		[DllImport("vocoder.dll")]
		static extern void mkvo(double[] gen, int allfrm, int MGCvecl, double[,] MGC, double[,] LF0, int BAPvecl, double[,] BAP, int stage, bool use_log_gain, int sampling_rate, int fperiod, double alpha, double beta, double[] volume);


		[DllImport("msvcrt.dll", CallingConvention = CallingConvention.Cdecl, SetLastError = false)]
		static extern int system(string command);


		static void pdwrt(int poz, string paty, HTSstm stm)
		{
   
			Directory.CreateDirectory(paty);
			using (FileStream fs = File.Create(paty + "\\dat.p")) {
				fs.Write(stm.pdbuf, poz, stm.pdsz);
			}
		}

		static void sinope(brnchy turi, string inpa, HTSstm stm)
		{
			turi.key.makefol(inpa);

			if (turi.pdn != -1) {
				pdwrt(turi.pdn, inpa + "\\-", stm);

			} else {
				sinope(turi.ifn, inpa + "\\-", stm);
			}


			if (turi.pdy != -1) {
				pdwrt(turi.pdy, inpa + "\\+", stm);

			} else {
				sinope(turi.ify, inpa + "\\+", stm);
			}


		}

		static void testdump(string bdir, HTSvoice homo)
		{
			for (int i = 0; i < 10; i++) {
				if (homo.stmz[i] != null) {
					int ly = homo.stmz[i].HTStr.Length;
					for (int j = 0; j < ly; j++) {
						sinope(homo.stmz[i].HTStr[j][0], bdir + homo.stmz[i].stmsig + "\\T" + j + '\\', homo.stmz[i]);
					}

				}
			}

		}
		
		
		static labio[][] mklabzz(int nlabz)
		{
			labio[][] labik = new labio[nlabz][];
			
			
			for (int i1 = 0; i1 < nlabz; i1++) {
				string[] plab = File.ReadAllLines("inin.txt");
				int spo = plab.Length;
				labio[] lab = new labio[spo];
			
				for (int i = 0; i < spo; i++) {
					lab[i] = new labio(plab[i]);
				}
			
				labik[i1] = lab;
			}
			return labik;
		}
		

		static void Main(string[] args)
		{

			rydprv();

			HTSvoice[] homo = new HTSvoice[1];
			
			homo[0] = new HTSvoice(@"Q:\b\8t1f0");
			//homo[1] = new HTSvoice(@"Q:\b\8t1f2");
			
			int nlabz = 1;
			labio[][] labik = mklabzz(nlabz);
			
			
			for (int i = 0; i < nlabz; i++) {
			
			
				sstream sstmj = new sstream(homo, labik[i]);
				sstmj.setspeed(1.0);
				//sstmj.mkdurseg(1.2);
				//sstmj.mkdurseg(0.1);
			

				double[][,] gstrm = sstmj.produce();
				
				//=====================
				/*
				for (int stt = 0; stt < 3; stt++) {
					for (int t = 0; t < gstrm[stt].GetLength(0); t++) {
						string siya = "stm" + stt + ", t" + t + ": " + '\t';
					
						for (int vv = 0; vv < gstrm[stt].GetLength(1); vv++) {
							siya += gstrm[stt][t, vv].ToString("#0.000000") + ", " + '\t';
						}
					
						Console.WriteLine(siya);
					
					}
				}
				*/
				//=====================
				
				double[,] MGC = gstrm[0];
				double[,] LF0 = gstrm[1];
				double[,] BAP = gstrm[2];
				int allfrm = MGC.GetLength(0);
				int MGCvecl = MGC.GetLength(1);
				int BAPvecl = BAP.GetLength(1);
				
				int stage = 0;
				bool use_log_gain = true;
				int sampling_rate = 48000;
				int fperiod = 240;
				
				
				double alpha = 0.55;
				double beta = 0;
				double[] gen = new double[allfrm * fperiod];
				double[] volume = new double[allfrm];
				for (int i1 = 0; i1 < allfrm; i1++) {
					volume[i1] = 1.0;
				}
				
				mkvo(gen, allfrm, MGCvecl, MGC, LF0, BAPvecl, BAP, stage, use_log_gain, sampling_rate, fperiod, alpha, beta, volume);
				
				short[] viv = new short[22 + allfrm * fperiod];
				int alyi = allfrm * fperiod;
				Parallel.For(0, alyi, i1 => {
					viv[22 + i1] = (short)gen[i1];
				});
				
				
				Console.WriteLine("omanko");
				
				
			}
			
                                                          
			return;
		}

	}
}