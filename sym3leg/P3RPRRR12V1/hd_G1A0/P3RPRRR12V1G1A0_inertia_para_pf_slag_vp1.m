% Calculate inertia matrix for parallel robot
% P3RPRRR12V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR12V1G1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:19
% EndTime: 2020-08-06 18:21:20
% DurationCPUTime: 0.55s
% Computational Cost: add. (1329->153), mult. (1446->266), div. (189->7), fcn. (945->18), ass. (0->107)
t372 = 2 * rSges(2,3);
t371 = 2 * m(3) * rSges(3,1) * rSges(3,2) - 2 * Icges(3,4);
t345 = (rSges(3,3) + pkin(5));
t333 = sin(qJ(3,3));
t339 = cos(qJ(3,3));
t301 = t339 * rSges(3,1) - t333 * rSges(3,2);
t370 = m(3) * t301;
t335 = sin(qJ(3,2));
t341 = cos(qJ(3,2));
t302 = t341 * rSges(3,1) - t335 * rSges(3,2);
t369 = m(3) * t302;
t337 = sin(qJ(3,1));
t343 = cos(qJ(3,1));
t303 = t343 * rSges(3,1) - t337 * rSges(3,2);
t368 = m(3) * t303;
t354 = 1 / pkin(3);
t367 = m(3) * t354;
t366 = (-pkin(1) - t345) * m(3);
t351 = rSges(3,2) ^ 2;
t353 = rSges(3,1) ^ 2;
t365 = ((t351 + t353) * m(3) + Icges(3,3)) * t354;
t311 = t333 * pkin(3) + qJ(2,3);
t308 = 0.1e1 / t311;
t326 = 0.1e1 / t333;
t364 = t308 * t326;
t312 = t335 * pkin(3) + qJ(2,2);
t309 = 0.1e1 / t312;
t327 = 0.1e1 / t335;
t363 = t309 * t327;
t313 = t337 * pkin(3) + qJ(2,1);
t310 = 0.1e1 / t313;
t328 = 0.1e1 / t337;
t362 = t310 * t328;
t361 = t326 * t339;
t360 = t327 * t341;
t359 = t328 * t343;
t358 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) + Icges(3,2) + Icges(1,3);
t355 = pkin(1) ^ 2;
t357 = t353 + t355 + (2 * pkin(1) + t345) * t345;
t356 = rSges(2,3) ^ 2 + t355 + (-2 * pkin(1) + rSges(2,2)) * rSges(2,2);
t349 = qJ(2,1) ^ 2;
t348 = qJ(2,2) ^ 2;
t347 = qJ(2,3) ^ 2;
t346 = m(2) + m(3);
t344 = cos(qJ(1,1));
t342 = cos(qJ(1,2));
t340 = cos(qJ(1,3));
t338 = sin(qJ(1,1));
t336 = sin(qJ(1,2));
t334 = sin(qJ(1,3));
t332 = legFrame(1,3);
t331 = legFrame(2,3);
t330 = legFrame(3,3);
t324 = pkin(1) + pkin(5) + pkin(6);
t322 = cos(t332);
t321 = cos(t331);
t320 = cos(t330);
t319 = sin(t332);
t318 = sin(t331);
t317 = sin(t330);
t306 = rSges(3,1) * t366 + Icges(3,5);
t305 = -rSges(3,2) * t366 - Icges(3,6);
t304 = (t351 - t353) * m(3) + Icges(3,1) - Icges(3,2);
t300 = t366 - (pkin(1) - rSges(2,2)) * m(2);
t299 = t319 * t344 + t322 * t338;
t298 = t318 * t342 + t321 * t336;
t297 = t317 * t340 + t320 * t334;
t296 = -t319 * t338 + t322 * t344;
t295 = -t318 * t336 + t321 * t342;
t294 = -t317 * t334 + t320 * t340;
t293 = -t313 * t344 + t324 * t338;
t292 = -t312 * t342 + t324 * t336;
t291 = -t311 * t340 + t324 * t334;
t290 = t338 * t313 + t324 * t344;
t289 = t336 * t312 + t324 * t342;
t288 = t334 * t311 + t324 * t340;
t287 = t337 * t305 + t306 * t343;
t286 = t335 * t305 + t306 * t341;
t285 = t333 * t305 + t306 * t339;
t284 = (-t303 * t367 + t343 * t346) * t328;
t283 = (-t302 * t367 + t341 * t346) * t327;
t282 = (-t301 * t367 + t339 * t346) * t326;
t281 = (-t287 * t354 + t300 * t343) * t328;
t280 = (-t286 * t354 + t300 * t341) * t327;
t279 = (-t285 * t354 + t300 * t339) * t326;
t278 = t319 * t290 + t293 * t322;
t277 = t318 * t289 + t292 * t321;
t276 = t317 * t288 + t291 * t320;
t275 = t290 * t322 - t319 * t293;
t274 = t289 * t321 - t318 * t292;
t273 = t288 * t320 - t317 * t291;
t272 = (qJ(2,1) * t372 + t349 + t356) * m(2) + (t304 * t343 + t337 * t371) * t343 + (t349 + 0.2e1 * (rSges(3,1) * t337 + rSges(3,2) * t343) * qJ(2,1) + t357) * m(3) + t358;
t271 = (qJ(2,2) * t372 + t348 + t356) * m(2) + (t304 * t341 + t335 * t371) * t341 + (t348 + 0.2e1 * (rSges(3,1) * t335 + rSges(3,2) * t341) * qJ(2,2) + t357) * m(3) + t358;
t270 = (qJ(2,3) * t372 + t347 + t356) * m(2) + (t304 * t339 + t333 * t371) * t339 + (t347 + 0.2e1 * (rSges(3,1) * t333 + rSges(3,2) * t339) * qJ(2,3) + t357) * m(3) + t358;
t269 = (t278 * t346 + t299 * t300) * t310;
t268 = (t277 * t346 + t298 * t300) * t309;
t267 = (t276 * t346 + t297 * t300) * t308;
t266 = (t275 * t346 + t296 * t300) * t310;
t265 = (t274 * t346 + t295 * t300) * t309;
t264 = (t273 * t346 + t294 * t300) * t308;
t263 = (t272 * t299 + t278 * t300) * t310;
t262 = (t271 * t298 + t277 * t300) * t309;
t261 = (t270 * t297 + t276 * t300) * t308;
t260 = (t272 * t296 + t275 * t300) * t310;
t259 = (t271 * t295 + t274 * t300) * t309;
t258 = (t270 * t294 + t273 * t300) * t308;
t1 = [m(4) + (t260 * t296 + t266 * t275) * t310 + (t259 * t295 + t265 * t274) * t309 + (t258 * t294 + t264 * t273) * t308, (t260 * t299 + t266 * t278) * t310 + (t259 * t298 + t265 * t277) * t309 + (t258 * t297 + t264 * t276) * t308, t264 * t361 + t265 * t360 + t266 * t359 + (-(t275 * t368 + t287 * t296) * t362 - (t274 * t369 + t286 * t295) * t363 - (t273 * t370 + t285 * t294) * t364) * t354; (t263 * t296 + t269 * t275) * t310 + (t262 * t295 + t268 * t274) * t309 + (t261 * t294 + t267 * t273) * t308, m(4) + (t263 * t299 + t269 * t278) * t310 + (t262 * t298 + t268 * t277) * t309 + (t261 * t297 + t267 * t276) * t308, t267 * t361 + t268 * t360 + t269 * t359 + (-(t278 * t368 + t287 * t299) * t362 - (t277 * t369 + t286 * t298) * t363 - (t276 * t370 + t285 * t297) * t364) * t354; (t275 * t284 + t281 * t296) * t310 + (t274 * t283 + t280 * t295) * t309 + (t273 * t282 + t279 * t294) * t308, (t278 * t284 + t281 * t299) * t310 + (t277 * t283 + t280 * t298) * t309 + (t276 * t282 + t279 * t297) * t308, m(4) + (t284 * t343 - (t343 * t368 - t365) * t354 * t328) * t328 + (t283 * t341 - (t341 * t369 - t365) * t354 * t327) * t327 + (t282 * t339 - (t339 * t370 - t365) * t354 * t326) * t326;];
MX  = t1;
