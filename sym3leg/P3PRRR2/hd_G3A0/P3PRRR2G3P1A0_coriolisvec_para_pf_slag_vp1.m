% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR2G3P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3P1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:01
% EndTime: 2020-03-09 21:20:02
% DurationCPUTime: 0.45s
% Computational Cost: add. (1828->97), mult. (1422->187), div. (1212->5), fcn. (1260->18), ass. (0->110)
t324 = cos(qJ(3,3));
t374 = t324 * pkin(1);
t325 = cos(qJ(3,2));
t373 = t325 * pkin(1);
t326 = cos(qJ(3,1));
t372 = t326 * pkin(1);
t328 = xDP(1);
t332 = 0.1e1 / pkin(2);
t334 = 0.1e1 / pkin(1);
t348 = t332 * t334;
t342 = t328 * t348;
t310 = -legFrame(3,2) + qJ(2,3);
t306 = qJ(3,3) + t310;
t300 = sin(t306);
t287 = pkin(2) * t300 + pkin(1) * sin(t310);
t321 = sin(qJ(3,3));
t315 = 0.1e1 / t321;
t368 = t287 * t315;
t278 = t342 * t368;
t327 = xDP(2);
t343 = t327 * t348;
t303 = cos(t306);
t290 = -pkin(2) * t303 - pkin(1) * cos(t310);
t365 = t290 * t315;
t281 = t343 * t365;
t269 = t278 + t281;
t349 = t328 * t334;
t350 = t327 * t334;
t356 = t303 * t315;
t359 = t300 * t315;
t275 = -t349 * t359 + t350 * t356;
t266 = t275 + t269;
t371 = t266 * t269;
t311 = -legFrame(2,2) + qJ(2,2);
t307 = qJ(3,2) + t311;
t301 = sin(t307);
t288 = pkin(2) * t301 + pkin(1) * sin(t311);
t322 = sin(qJ(3,2));
t316 = 0.1e1 / t322;
t367 = t288 * t316;
t279 = t342 * t367;
t304 = cos(t307);
t291 = -pkin(2) * t304 - pkin(1) * cos(t311);
t364 = t291 * t316;
t282 = t343 * t364;
t270 = t279 + t282;
t355 = t304 * t316;
t358 = t301 * t316;
t276 = -t349 * t358 + t350 * t355;
t267 = t276 + t270;
t370 = t267 * t270;
t312 = -legFrame(1,2) + qJ(2,1);
t308 = qJ(3,1) + t312;
t302 = sin(t308);
t289 = pkin(2) * t302 + pkin(1) * sin(t312);
t323 = sin(qJ(3,1));
t317 = 0.1e1 / t323;
t366 = t289 * t317;
t280 = t342 * t366;
t305 = cos(t308);
t292 = -pkin(2) * t305 - pkin(1) * cos(t312);
t363 = t292 * t317;
t283 = t343 * t363;
t271 = t280 + t283;
t354 = t305 * t317;
t357 = t302 * t317;
t277 = -t349 * t357 + t350 * t354;
t268 = t277 + t271;
t369 = t268 * t271;
t362 = (t321 * rSges(3,1) + t324 * rSges(3,2)) * t315;
t361 = (t322 * rSges(3,1) + t325 * rSges(3,2)) * t316;
t360 = (t323 * rSges(3,1) + t326 * rSges(3,2)) * t317;
t353 = t315 * t334;
t352 = t316 * t334;
t351 = t317 * t334;
t347 = 0.2e1 * pkin(1) * pkin(2);
t313 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t346 = t275 ^ 2 * t362;
t345 = t276 ^ 2 * t361;
t344 = t277 ^ 2 * t360;
t263 = t278 / 0.2e1 + t281 / 0.2e1 + t275;
t341 = t263 * t269 * t362;
t264 = t279 / 0.2e1 + t282 / 0.2e1 + t276;
t340 = t264 * t270 * t361;
t265 = t280 / 0.2e1 + t283 / 0.2e1 + t277;
t339 = t265 * t271 * t360;
t338 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(2,3) - Icges(3,3);
t337 = (-rSges(3,1) * t324 + rSges(3,2) * t321) * pkin(1);
t336 = (-rSges(3,1) * t325 + rSges(3,2) * t322) * pkin(1);
t335 = (-rSges(3,1) * t326 + rSges(3,2) * t323) * pkin(1);
t333 = pkin(1) ^ 2;
t331 = pkin(2) ^ 2;
t309 = t333 + t313;
t299 = -t313 * m(3) - Icges(3,3);
t286 = -Icges(3,3) + (-t313 + t335) * m(3);
t285 = -Icges(3,3) + (-t313 + t336) * m(3);
t284 = -Icges(3,3) + (-t313 + t337) * m(3);
t262 = ((-pkin(2) * t268 - t277 * t372) * t277 - pkin(2) * t369) * t351;
t261 = ((-pkin(2) * t267 - t276 * t373) * t276 - pkin(2) * t370) * t352;
t260 = ((-pkin(2) * t266 - t275 * t374) * t275 - pkin(2) * t371) * t353;
t259 = ((t265 * t326 * t347 + t268 * t331 + t333 * t277) * t332 * t277 + (pkin(2) + t372) * t369) * t351;
t258 = ((t264 * t325 * t347 + t267 * t331 + t333 * t276) * t332 * t276 + (pkin(2) + t373) * t370) * t352;
t257 = ((t263 * t324 * t347 + t266 * t331 + t333 * t275) * t332 * t275 + (pkin(2) + t374) * t371) * t353;
t256 = t299 * t259 + t286 * t262;
t255 = t299 * t258 + t285 * t261;
t254 = t299 * t257 + t284 * t260;
t253 = t286 * t259 + (t338 + (-t309 + 0.2e1 * t335) * m(3)) * t262;
t252 = t285 * t258 + (t338 + (-t309 + 0.2e1 * t336) * m(3)) * t261;
t251 = t284 * t257 + (t338 + (-t309 + 0.2e1 * t337) * m(3)) * t260;
t1 = [(-t251 * t359 - t252 * t358 - t253 * t357 + (t254 * t368 + t255 * t367 + t256 * t366) * t332) * t334 + (0.2e1 * t300 * t341 + 0.2e1 * t301 * t340 + 0.2e1 * t302 * t339 + (t287 * t346 + t288 * t345 + t289 * t344) * t332) * m(3); (t251 * t356 + t252 * t355 + t253 * t354 + (t254 * t365 + t255 * t364 + t256 * t363) * t332) * t334 + (-0.2e1 * t303 * t341 - 0.2e1 * t304 * t340 - 0.2e1 * t305 * t339 + (t290 * t346 + t291 * t345 + t292 * t344) * t332) * m(3); 0;];
taucX  = t1;
