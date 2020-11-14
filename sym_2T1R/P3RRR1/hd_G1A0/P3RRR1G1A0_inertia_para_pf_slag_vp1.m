% Calculate inertia matrix for parallel robot
% P3RRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 15:38
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRR1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRR1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RRR1G1A0_inertia_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRR1G1A0_inertia_para_pf_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RRR1G1A0_inertia_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3RRR1G1A0_inertia_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3RRR1G1A0_inertia_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRR1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRR1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 15:38:14
% EndTime: 2019-05-03 15:38:15
% DurationCPUTime: 0.53s
% Computational Cost: add. (2596->132), mult. (3201->254), div. (540->5), fcn. (2870->20), ass. (0->118)
t321 = qJ(1,3) + qJ(2,3);
t307 = sin(t321);
t310 = cos(t321);
t327 = sin(qJ(1,3));
t330 = cos(qJ(1,3));
t374 = 0.1e1 / (t307 * t330 - t327 * t310);
t322 = qJ(1,2) + qJ(2,2);
t308 = sin(t322);
t311 = cos(t322);
t328 = sin(qJ(1,2));
t331 = cos(qJ(1,2));
t373 = 0.1e1 / (t308 * t331 - t328 * t311);
t323 = qJ(1,1) + qJ(2,1);
t309 = sin(t323);
t312 = cos(t323);
t329 = sin(qJ(1,1));
t332 = cos(qJ(1,1));
t372 = 0.1e1 / (t309 * t332 - t329 * t312);
t371 = m(2) * pkin(1);
t324 = legFrame(3,3);
t313 = sin(t324);
t316 = cos(t324);
t274 = t316 * t307 + t313 * t310;
t268 = pkin(1) * (t313 * t330 + t316 * t327) + t274 * pkin(2);
t370 = t268 * t374;
t325 = legFrame(2,3);
t314 = sin(t325);
t317 = cos(t325);
t276 = t317 * t308 + t314 * t311;
t269 = pkin(1) * (t314 * t331 + t317 * t328) + t276 * pkin(2);
t369 = t269 * t373;
t326 = legFrame(1,3);
t315 = sin(t326);
t318 = cos(t326);
t278 = t318 * t309 + t315 * t312;
t270 = pkin(1) * (t315 * t332 + t318 * t329) + t278 * pkin(2);
t368 = t270 * t372;
t275 = -t307 * t313 + t316 * t310;
t271 = pkin(1) * (-t327 * t313 + t316 * t330) + t275 * pkin(2);
t367 = t271 * t374;
t277 = -t308 * t314 + t317 * t311;
t272 = pkin(1) * (-t328 * t314 + t317 * t331) + t277 * pkin(2);
t366 = t272 * t373;
t279 = -t309 * t315 + t318 * t312;
t273 = pkin(1) * (-t329 * t315 + t318 * t332) + t279 * pkin(2);
t365 = t273 * t372;
t364 = t274 * t374;
t363 = t275 * t374;
t362 = t276 * t373;
t361 = t277 * t373;
t360 = t278 * t372;
t359 = t279 * t372;
t356 = rSges(2,1) ^ 2 + rSges(2,2) ^ 2;
t303 = t356 * m(2) + Icges(2,3);
t344 = 0.1e1 / pkin(2);
t358 = t303 * t344;
t345 = 0.1e1 / pkin(1);
t357 = t344 * t345;
t348 = (-(-rSges(2,1) * t327 + rSges(2,2) * t330) * t307 + (rSges(2,1) * t330 + rSges(2,2) * t327) * t310) * t371;
t265 = t348 + t303;
t355 = t265 * t374 * t344;
t347 = (-(-rSges(2,1) * t328 + rSges(2,2) * t331) * t308 + (rSges(2,1) * t331 + rSges(2,2) * t328) * t311) * t371;
t266 = t347 + t303;
t354 = t266 * t373 * t344;
t346 = (-(-rSges(2,1) * t329 + rSges(2,2) * t332) * t309 + (rSges(2,1) * t332 + rSges(2,2) * t329) * t312) * t371;
t267 = t346 + t303;
t353 = t267 * t372 * t344;
t352 = t374 * t358;
t351 = t373 * t358;
t350 = t372 * t358;
t349 = Icges(1,3) + Icges(2,3) + (pkin(1) ^ 2 + t356) * m(2) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1);
t341 = koppelP(1,1);
t340 = koppelP(2,1);
t339 = koppelP(3,1);
t338 = koppelP(1,2);
t337 = koppelP(2,2);
t336 = koppelP(3,2);
t335 = rSges(3,1);
t334 = rSges(3,2);
t333 = xP(3);
t320 = cos(t333);
t319 = sin(t333);
t296 = -t319 * t338 + t320 * t341;
t295 = -t319 * t337 + t320 * t340;
t294 = -t319 * t336 + t320 * t339;
t293 = -t319 * t341 - t320 * t338;
t292 = -t319 * t340 - t320 * t337;
t291 = -t319 * t339 - t320 * t336;
t290 = m(3) * (-t319 * t334 + t320 * t335);
t289 = m(3) * (-t319 * t335 - t320 * t334);
t264 = 0.2e1 * t346 + t349;
t263 = 0.2e1 * t347 + t349;
t262 = 0.2e1 * t348 + t349;
t261 = (t278 * t296 + t279 * t293) * t345 * t372;
t260 = (t276 * t295 + t277 * t292) * t345 * t373;
t259 = (t274 * t294 + t275 * t291) * t345 * t374;
t258 = (t270 * t296 + t273 * t293) * t372 * t357;
t257 = (t269 * t295 + t272 * t292) * t373 * t357;
t256 = (t268 * t294 + t271 * t291) * t374 * t357;
t255 = (t267 * t359 - t273 * t350) * t345;
t254 = (t267 * t360 - t270 * t350) * t345;
t253 = (t266 * t361 - t272 * t351) * t345;
t252 = (t266 * t362 - t269 * t351) * t345;
t251 = (t265 * t363 - t271 * t352) * t345;
t250 = (t265 * t364 - t268 * t352) * t345;
t249 = (t264 * t359 - t273 * t353) * t345;
t248 = (t264 * t360 - t270 * t353) * t345;
t247 = (t263 * t361 - t272 * t354) * t345;
t246 = (t263 * t362 - t269 * t354) * t345;
t245 = (t262 * t363 - t271 * t355) * t345;
t244 = (t262 * t364 - t268 * t355) * t345;
t243 = -t258 * t303 + t261 * t267;
t242 = -t257 * t303 + t260 * t266;
t241 = -t256 * t303 + t259 * t265;
t240 = -t258 * t267 + t261 * t264;
t239 = -t257 * t266 + t260 * t263;
t238 = -t256 * t265 + t259 * t262;
t1 = [m(3) + (t245 * t363 + t247 * t361 + t249 * t359 + (-t251 * t367 - t253 * t366 - t255 * t365) * t344) * t345, (t245 * t364 + t247 * t362 + t249 * t360 + (-t251 * t370 - t253 * t369 - t255 * t368) * t344) * t345, t245 * t259 + t247 * t260 + t249 * t261 - t251 * t256 - t253 * t257 - t255 * t258 + t289; (t244 * t363 + t246 * t361 + t248 * t359 + (-t250 * t367 - t252 * t366 - t254 * t365) * t344) * t345, m(3) + (t244 * t364 + t246 * t362 + t248 * t360 + (-t250 * t370 - t252 * t369 - t254 * t368) * t344) * t345, t244 * t259 + t246 * t260 + t248 * t261 - t250 * t256 - t252 * t257 - t254 * t258 + t290; t289 + (t238 * t363 + t239 * t361 + t240 * t359 + (-t241 * t367 - t242 * t366 - t243 * t365) * t344) * t345, t290 + (t238 * t364 + t239 * t362 + t240 * t360 + (-t241 * t370 - t242 * t369 - t243 * t368) * t344) * t345, t240 * t261 - t243 * t258 + t239 * t260 - t242 * t257 + t238 * t259 - t241 * t256 + Icges(3,3) + m(3) * (t334 ^ 2 + t335 ^ 2);];
MX  = t1;
