% Calculate inertia matrix for parallel robot
% P3PRRR1G3P3A0
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
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR1G3P3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:45
% EndTime: 2020-03-09 21:06:45
% DurationCPUTime: 0.56s
% Computational Cost: add. (3054->114), mult. (2076->226), div. (405->5), fcn. (1503->18), ass. (0->115)
t299 = pkin(7) + qJ(2,3);
t290 = qJ(3,3) + t299;
t276 = sin(t290);
t279 = cos(t290);
t284 = sin(t299);
t287 = cos(t299);
t362 = 0.1e1 / (t276 * t287 - t284 * t279);
t300 = pkin(7) + qJ(2,2);
t291 = qJ(3,2) + t300;
t277 = sin(t291);
t280 = cos(t291);
t285 = sin(t300);
t288 = cos(t300);
t361 = 0.1e1 / (t277 * t288 - t285 * t280);
t301 = pkin(7) + qJ(2,1);
t292 = qJ(3,1) + t301;
t278 = sin(t292);
t281 = cos(t292);
t286 = sin(t301);
t289 = cos(t301);
t360 = 0.1e1 / (t278 * t289 - t286 * t281);
t359 = m(3) * pkin(2);
t336 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t274 = t336 * m(3) + Icges(3,3);
t312 = ((t284 * rSges(3,1) - rSges(3,2) * t287) * t276 + (rSges(3,1) * t287 + t284 * rSges(3,2)) * t279) * t359;
t247 = t312 + t274;
t308 = 0.1e1 / pkin(3);
t358 = t247 * t308;
t311 = ((t285 * rSges(3,1) - rSges(3,2) * t288) * t277 + (rSges(3,1) * t288 + t285 * rSges(3,2)) * t280) * t359;
t248 = t311 + t274;
t357 = t248 * t308;
t310 = ((t286 * rSges(3,1) - rSges(3,2) * t289) * t278 + (rSges(3,1) * t289 + t286 * rSges(3,2)) * t281) * t359;
t249 = t310 + t274;
t356 = t249 * t308;
t259 = pkin(2) * t284 + pkin(3) * t276;
t355 = t362 * t259;
t354 = t362 * t276;
t353 = t362 * t279;
t303 = legFrame(3,2);
t293 = sin(t303);
t352 = t362 * t293;
t309 = 0.1e1 / pkin(2);
t351 = t362 * t309;
t260 = pkin(2) * t285 + pkin(3) * t277;
t350 = t361 * t260;
t349 = t361 * t277;
t348 = t361 * t280;
t304 = legFrame(2,2);
t294 = sin(t304);
t347 = t361 * t294;
t346 = t361 * t309;
t261 = pkin(2) * t286 + pkin(3) * t278;
t345 = t360 * t261;
t344 = t360 * t278;
t343 = t360 * t281;
t305 = legFrame(1,2);
t295 = sin(t305);
t342 = t360 * t295;
t341 = t360 * t309;
t340 = t274 * t308;
t296 = cos(t303);
t339 = t296 * t309;
t297 = cos(t304);
t338 = t297 * t309;
t298 = cos(t305);
t337 = t298 * t309;
t262 = pkin(2) * t287 + pkin(3) * t279;
t335 = t262 * t358;
t263 = pkin(2) * t288 + pkin(3) * t280;
t334 = t263 * t357;
t264 = pkin(2) * t289 + pkin(3) * t281;
t333 = t264 * t356;
t332 = t262 * t352;
t331 = t279 * t352;
t330 = t296 * t353;
t329 = t293 * t351;
t328 = t263 * t347;
t327 = t280 * t347;
t326 = t297 * t348;
t325 = t294 * t346;
t324 = t264 * t342;
t323 = t281 * t342;
t322 = t298 * t343;
t321 = t295 * t341;
t320 = t362 * t262 * t296;
t319 = t361 * t263 * t297;
t318 = t360 * t264 * t298;
t317 = t262 * t340;
t316 = t263 * t340;
t315 = t264 * t340;
t302 = m(1) + m(2) + m(3);
t314 = (t293 * t296 + t294 * t297 + t295 * t298) * t302;
t313 = Icges(2,3) + Icges(3,3) + (pkin(2) ^ 2 + t336) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t246 = 0.2e1 * t310 + t313;
t245 = 0.2e1 * t311 + t313;
t244 = 0.2e1 * t312 + t313;
t243 = (-t249 * t278 + t261 * t340) * t341;
t242 = (-t248 * t277 + t260 * t340) * t346;
t241 = (-t247 * t276 + t259 * t340) * t351;
t240 = (t249 * t343 - t315 * t360) * t337;
t239 = (t248 * t348 - t316 * t361) * t338;
t238 = (t247 * t353 - t317 * t362) * t339;
t237 = (-t249 * t281 + t315) * t321;
t236 = (-t248 * t280 + t316) * t325;
t235 = (-t247 * t279 + t317) * t329;
t234 = (-t246 * t278 + t261 * t356) * t341;
t233 = (-t245 * t277 + t260 * t357) * t346;
t232 = (-t244 * t276 + t259 * t358) * t351;
t231 = (t246 * t343 - t333 * t360) * t337;
t230 = (t245 * t348 - t334 * t361) * t338;
t229 = (t244 * t353 - t335 * t362) * t339;
t228 = (-t246 * t281 + t333) * t321;
t227 = (-t245 * t280 + t334) * t325;
t226 = (-t244 * t279 + t335) * t329;
t1 = [m(4) + (t293 ^ 2 + t294 ^ 2 + t295 ^ 2) * t302 + (t229 * t330 + t230 * t326 + t231 * t322 + (-t238 * t320 - t239 * t319 - t240 * t318) * t308) * t309, (-t229 * t331 - t230 * t327 - t231 * t323 + (t238 * t332 + t239 * t328 + t240 * t324) * t308) * t309 + t314, (-t229 * t354 - t230 * t349 - t231 * t344 + (t238 * t355 + t239 * t350 + t240 * t345) * t308) * t309; (t226 * t330 + t227 * t326 + t228 * t322 + (-t235 * t320 - t236 * t319 - t237 * t318) * t308) * t309 + t314, m(4) + (t296 ^ 2 + t297 ^ 2 + t298 ^ 2) * t302 + (-t226 * t331 - t227 * t327 - t228 * t323 + (t235 * t332 + t236 * t328 + t237 * t324) * t308) * t309, (-t226 * t354 - t227 * t349 - t228 * t344 + (t235 * t355 + t236 * t350 + t237 * t345) * t308) * t309; (t232 * t330 + t233 * t326 + t234 * t322 + (-t241 * t320 - t242 * t319 - t243 * t318) * t308) * t309, (-t232 * t331 - t233 * t327 - t234 * t323 + (t241 * t332 + t242 * t328 + t243 * t324) * t308) * t309, m(4) + (-t232 * t354 - t233 * t349 - t234 * t344 + (t241 * t355 + t242 * t350 + t243 * t345) * t308) * t309;];
MX  = t1;
