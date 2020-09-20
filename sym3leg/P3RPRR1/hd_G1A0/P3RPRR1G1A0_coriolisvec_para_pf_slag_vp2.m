% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:23
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G1P1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:23:14
% EndTime: 2020-03-09 21:23:15
% DurationCPUTime: 0.68s
% Computational Cost: add. (4756->137), mult. (3372->218), div. (840->7), fcn. (2532->32), ass. (0->100)
t361 = 2 * pkin(2);
t332 = sin(pkin(7));
t371 = pkin(1) * t332;
t333 = cos(pkin(7));
t370 = pkin(2) * t333;
t369 = m(3) * pkin(2) + mrSges(2,1);
t325 = legFrame(3,3) + qJ(1,3);
t318 = pkin(7) + t325;
t315 = qJ(3,3) + t318;
t309 = sin(t315);
t292 = -pkin(1) * sin(t325) - pkin(2) * sin(t318) - pkin(3) * t309;
t329 = pkin(7) + qJ(3,3);
t337 = sin(qJ(3,3));
t306 = 0.1e1 / (pkin(1) * sin(t329) + t337 * pkin(2));
t343 = xDP(2);
t348 = 0.1e1 / pkin(3);
t365 = t343 * t348;
t280 = t292 * t306 * t365;
t312 = cos(t315);
t295 = -pkin(1) * cos(t325) - pkin(2) * cos(t318) - pkin(3) * t312;
t344 = xDP(1);
t364 = t344 * t348;
t283 = t295 * t306 * t364;
t271 = t283 + t280;
t277 = (t309 * t343 + t312 * t344) * t306;
t268 = t277 + t271;
t368 = t268 * t271;
t326 = legFrame(2,3) + qJ(1,2);
t319 = pkin(7) + t326;
t316 = qJ(3,2) + t319;
t310 = sin(t316);
t293 = -pkin(1) * sin(t326) - pkin(2) * sin(t319) - pkin(3) * t310;
t330 = pkin(7) + qJ(3,2);
t338 = sin(qJ(3,2));
t307 = 0.1e1 / (pkin(1) * sin(t330) + t338 * pkin(2));
t281 = t293 * t307 * t365;
t313 = cos(t316);
t296 = -pkin(1) * cos(t326) - pkin(2) * cos(t319) - pkin(3) * t313;
t284 = t296 * t307 * t364;
t272 = t284 + t281;
t278 = (t310 * t343 + t313 * t344) * t307;
t269 = t278 + t272;
t367 = t269 * t272;
t327 = legFrame(1,3) + qJ(1,1);
t320 = pkin(7) + t327;
t317 = qJ(3,1) + t320;
t311 = sin(t317);
t294 = -pkin(1) * sin(t327) - pkin(2) * sin(t320) - pkin(3) * t311;
t331 = pkin(7) + qJ(3,1);
t339 = sin(qJ(3,1));
t308 = 0.1e1 / (pkin(1) * sin(t331) + t339 * pkin(2));
t282 = t294 * t308 * t365;
t314 = cos(t317);
t297 = -pkin(1) * cos(t327) - pkin(2) * cos(t320) - pkin(3) * t314;
t285 = t297 * t308 * t364;
t273 = t285 + t282;
t279 = (t311 * t343 + t314 * t344) * t308;
t270 = t279 + t273;
t366 = t270 * t273;
t363 = pkin(3) * t361;
t362 = 0.2e1 * pkin(1);
t321 = t333 * pkin(1) + pkin(2);
t340 = cos(qJ(3,3));
t360 = -t340 * mrSges(3,1) + mrSges(3,2) * t337;
t341 = cos(qJ(3,2));
t359 = -t341 * mrSges(3,1) + mrSges(3,2) * t338;
t342 = cos(qJ(3,1));
t358 = -t342 * mrSges(3,1) + mrSges(3,2) * t339;
t265 = t283 / 0.2e1 + t280 / 0.2e1 + t277;
t322 = cos(t329);
t349 = pkin(2) ^ 2;
t350 = pkin(1) ^ 2;
t328 = t349 + t350;
t347 = pkin(3) ^ 2;
t259 = (t265 * t340 * t363 + t268 * t347 + t277 * t328 + (pkin(3) * t265 * t322 + t277 * t370) * t362) * t348 * t306 * t277 + (t321 * t340 - t337 * t371 + pkin(3)) / (t321 * t337 + t340 * t371) * t368;
t262 = (-pkin(3) * t368 + (-pkin(3) * t268 + (-pkin(1) * t322 - t340 * pkin(2)) * t277) * t277) * t306;
t304 = mrSges(3,1) * t371 + t321 * mrSges(3,2);
t305 = t321 * mrSges(3,1) - mrSges(3,2) * t371;
t286 = t304 * t337 - t305 * t340 - Ifges(3,3);
t289 = t304 * t340 + t337 * t305;
t357 = (t277 ^ 2 * t289 - Ifges(3,3) * t259 + t286 * t262) * t306;
t266 = t284 / 0.2e1 + t281 / 0.2e1 + t278;
t323 = cos(t330);
t260 = (t266 * t341 * t363 + t269 * t347 + t278 * t328 + (pkin(3) * t266 * t323 + t278 * t370) * t362) * t348 * t307 * t278 + (t321 * t341 - t338 * t371 + pkin(3)) / (t321 * t338 + t341 * t371) * t367;
t263 = (-pkin(3) * t367 + (-pkin(3) * t269 + (-pkin(1) * t323 - t341 * pkin(2)) * t278) * t278) * t307;
t287 = t304 * t338 - t305 * t341 - Ifges(3,3);
t290 = t304 * t341 + t338 * t305;
t356 = (t278 ^ 2 * t290 - Ifges(3,3) * t260 + t287 * t263) * t307;
t267 = t285 / 0.2e1 + t282 / 0.2e1 + t279;
t324 = cos(t331);
t261 = (t267 * t342 * t363 + t270 * t347 + t279 * t328 + (pkin(3) * t267 * t324 + t279 * t370) * t362) * t348 * t308 * t279 + (t342 * t321 - t339 * t371 + pkin(3)) / (t339 * t321 + t342 * t371) * t366;
t264 = (-pkin(3) * t366 + (-t270 * pkin(3) + (-pkin(1) * t324 - t342 * pkin(2)) * t279) * t279) * t308;
t288 = t304 * t339 - t305 * t342 - Ifges(3,3);
t291 = t304 * t342 + t339 * t305;
t355 = (t279 ^ 2 * t291 - Ifges(3,3) * t261 + t288 * t264) * t308;
t351 = -(m(3) * t349) - (m(2) + m(3)) * t350 - Ifges(1,3) - Ifges(2,3) - Ifges(3,3);
t354 = t306 * (-0.2e1 * t265 * t271 * t289 + (t360 * t361 + (-(-t360 + t369) * t333 + (mrSges(3,1) * t337 + t340 * mrSges(3,2) + mrSges(2,2)) * t332) * t362 + t351) * t262 + t286 * t259);
t353 = t307 * (-0.2e1 * t266 * t272 * t290 + (t359 * t361 + (-(-t359 + t369) * t333 + (mrSges(3,1) * t338 + t341 * mrSges(3,2) + mrSges(2,2)) * t332) * t362 + t351) * t263 + t287 * t260);
t352 = t308 * (-0.2e1 * t267 * t273 * t291 + (t358 * t361 + (-(-t358 + t369) * t333 + (mrSges(3,1) * t339 + t342 * mrSges(3,2) + mrSges(2,2)) * t332) * t362 + t351) * t264 + t288 * t261);
t1 = [t314 * t352 + t313 * t353 + t312 * t354 + (t295 * t357 + t296 * t356 + t297 * t355) * t348; t311 * t352 + t310 * t353 + t309 * t354 + (t292 * t357 + t293 * t356 + t294 * t355) * t348; 0;];
taucX  = t1;
