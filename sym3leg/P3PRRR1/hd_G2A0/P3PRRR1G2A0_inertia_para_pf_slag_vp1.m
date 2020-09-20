% Calculate inertia matrix for parallel robot
% P3PRRR1G2P2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR1G2P2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2P2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:17
% EndTime: 2020-03-09 21:18:17
% DurationCPUTime: 0.56s
% Computational Cost: add. (3054->114), mult. (2076->223), div. (405->5), fcn. (1503->18), ass. (0->109)
t312 = pkin(7) + qJ(2,3);
t303 = qJ(3,3) + t312;
t289 = sin(t303);
t292 = cos(t303);
t297 = sin(t312);
t300 = cos(t312);
t372 = 0.1e1 / (-t300 * t289 + t292 * t297);
t313 = pkin(7) + qJ(2,2);
t304 = qJ(3,2) + t313;
t290 = sin(t304);
t293 = cos(t304);
t298 = sin(t313);
t301 = cos(t313);
t371 = 0.1e1 / (-t301 * t290 + t293 * t298);
t314 = pkin(7) + qJ(2,1);
t305 = qJ(3,1) + t314;
t291 = sin(t305);
t294 = cos(t305);
t299 = sin(t314);
t302 = cos(t314);
t370 = 0.1e1 / (-t302 * t291 + t294 * t299);
t369 = m(3) * pkin(2);
t368 = t372 * (pkin(2) * t300 + pkin(3) * t292);
t316 = legFrame(3,2);
t306 = sin(t316);
t367 = t372 * t306;
t366 = t371 * (pkin(2) * t301 + pkin(3) * t293);
t317 = legFrame(2,2);
t307 = sin(t317);
t365 = t371 * t307;
t364 = t370 * (pkin(2) * t302 + pkin(3) * t294);
t318 = legFrame(1,2);
t308 = sin(t318);
t363 = t370 * t308;
t362 = t372 * t289;
t361 = t372 * t292;
t360 = t371 * t290;
t359 = t371 * t293;
t358 = t370 * t291;
t357 = t370 * t294;
t352 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t287 = t352 * m(3) + Icges(3,3);
t321 = 0.1e1 / pkin(3);
t356 = t287 * t321;
t309 = cos(t316);
t322 = 0.1e1 / pkin(2);
t355 = t309 * t322;
t310 = cos(t317);
t354 = t310 * t322;
t311 = cos(t318);
t353 = t311 * t322;
t325 = ((t297 * rSges(3,1) - rSges(3,2) * t300) * t289 + (rSges(3,1) * t300 + t297 * rSges(3,2)) * t292) * t369;
t260 = t325 + t287;
t272 = pkin(2) * t297 + pkin(3) * t289;
t351 = t260 * t272 * t321;
t324 = ((t298 * rSges(3,1) - rSges(3,2) * t301) * t290 + (rSges(3,1) * t301 + t298 * rSges(3,2)) * t293) * t369;
t261 = t324 + t287;
t273 = pkin(2) * t298 + pkin(3) * t290;
t350 = t261 * t273 * t321;
t323 = ((t299 * rSges(3,1) - rSges(3,2) * t302) * t291 + (rSges(3,1) * t302 + t299 * rSges(3,2)) * t294) * t369;
t262 = t323 + t287;
t274 = pkin(2) * t299 + pkin(3) * t291;
t349 = t262 * t274 * t321;
t348 = t272 * t367;
t347 = t372 * t272 * t309;
t346 = t321 * t368;
t345 = t289 * t367;
t344 = t322 * t367;
t343 = t273 * t365;
t342 = t371 * t273 * t310;
t341 = t321 * t366;
t340 = t290 * t365;
t339 = t322 * t365;
t338 = t274 * t363;
t337 = t370 * t274 * t311;
t336 = t321 * t364;
t335 = t291 * t363;
t334 = t322 * t363;
t333 = t309 * t362;
t332 = t310 * t360;
t331 = t311 * t358;
t330 = t272 * t356;
t329 = t273 * t356;
t328 = t274 * t356;
t315 = m(1) + m(2) + m(3);
t327 = (t306 * t309 + t307 * t310 + t308 * t311) * t315;
t326 = Icges(2,3) + Icges(3,3) + (pkin(2) ^ 2 + t352) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t259 = 0.2e1 * t323 + t326;
t258 = 0.2e1 * t324 + t326;
t257 = 0.2e1 * t325 + t326;
t256 = (-t262 * t357 + t287 * t336) * t322;
t255 = (-t261 * t359 + t287 * t341) * t322;
t254 = (-t260 * t361 + t287 * t346) * t322;
t253 = (t262 * t291 - t328) * t334;
t252 = (t261 * t290 - t329) * t339;
t251 = (t260 * t289 - t330) * t344;
t250 = (-t262 * t358 + t328 * t370) * t353;
t249 = (-t261 * t360 + t329 * t371) * t354;
t248 = (-t260 * t362 + t330 * t372) * t355;
t247 = (-t259 * t357 + t262 * t336) * t322;
t246 = (-t258 * t359 + t261 * t341) * t322;
t245 = (-t257 * t361 + t260 * t346) * t322;
t244 = (t259 * t291 - t349) * t334;
t243 = (t258 * t290 - t350) * t339;
t242 = (t257 * t289 - t351) * t344;
t241 = (-t259 * t358 + t349 * t370) * t353;
t240 = (-t258 * t360 + t350 * t371) * t354;
t239 = (-t257 * t362 + t351 * t372) * t355;
t1 = [m(4) + (t306 ^ 2 + t307 ^ 2 + t308 ^ 2) * t315 + (-t239 * t333 - t240 * t332 - t241 * t331 + (t248 * t347 + t249 * t342 + t250 * t337) * t321) * t322, (t239 * t345 + t240 * t340 + t241 * t335 + (-t248 * t348 - t249 * t343 - t250 * t338) * t321) * t322 + t327, (-t239 * t361 - t240 * t359 - t241 * t357 + (t248 * t368 + t249 * t366 + t250 * t364) * t321) * t322; (-t242 * t333 - t243 * t332 - t244 * t331 + (t251 * t347 + t252 * t342 + t253 * t337) * t321) * t322 + t327, m(4) + (t309 ^ 2 + t310 ^ 2 + t311 ^ 2) * t315 + (t242 * t345 + t243 * t340 + t244 * t335 + (-t251 * t348 - t252 * t343 - t253 * t338) * t321) * t322, (-t242 * t361 - t243 * t359 - t244 * t357 + (t251 * t368 + t252 * t366 + t253 * t364) * t321) * t322; (-t245 * t333 - t246 * t332 - t247 * t331 + (t254 * t347 + t255 * t342 + t256 * t337) * t321) * t322, (t245 * t345 + t246 * t340 + t247 * t335 + (-t254 * t348 - t255 * t343 - t256 * t338) * t321) * t322, m(4) + (-t245 * t361 - t246 * t359 - t247 * t357 + (t254 * t368 + t255 * t366 + t256 * t364) * t321) * t322;];
MX  = t1;
