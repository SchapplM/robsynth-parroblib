% Calculate inertia matrix for parallel robot
% P3RPRR1G3A0
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRR1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRR1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:26:51
% EndTime: 2020-03-09 21:26:52
% DurationCPUTime: 0.77s
% Computational Cost: add. (2637->178), mult. (2346->234), div. (243->4), fcn. (1422->62), ass. (0->130)
t341 = sin(qJ(3,3));
t358 = pkin(7) + qJ(3,3);
t290 = 0.1e1 / (pkin(1) * sin(t358) + t341 * pkin(2));
t370 = t290 / 0.2e1;
t342 = sin(qJ(3,2));
t359 = pkin(7) + qJ(3,2);
t291 = 0.1e1 / (pkin(1) * sin(t359) + t342 * pkin(2));
t368 = t291 / 0.2e1;
t343 = sin(qJ(3,1));
t360 = pkin(7) + qJ(3,1);
t292 = 0.1e1 / (pkin(1) * sin(t360) + t343 * pkin(2));
t366 = t292 / 0.2e1;
t351 = 0.1e1 / pkin(3);
t402 = t351 / 0.2e1;
t340 = legFrame(1,2);
t363 = qJ(1,1) + pkin(7);
t318 = t340 + t363;
t312 = qJ(3,1) + t318;
t319 = -t340 + t363;
t313 = qJ(3,1) + t319;
t287 = cos(t313) + cos(t312);
t401 = t287 * t366;
t339 = legFrame(2,2);
t362 = qJ(1,2) + pkin(7);
t316 = t339 + t362;
t310 = qJ(3,2) + t316;
t317 = -t339 + t362;
t311 = qJ(3,2) + t317;
t286 = cos(t311) + cos(t310);
t400 = t286 * t368;
t338 = legFrame(3,2);
t361 = qJ(1,3) + pkin(7);
t314 = t338 + t361;
t308 = qJ(3,3) + t314;
t315 = -t338 + t361;
t309 = qJ(3,3) + t315;
t285 = cos(t309) + cos(t308);
t399 = t285 * t370;
t284 = -sin(t312) + sin(t313);
t398 = t284 * t366;
t283 = -sin(t310) + sin(t311);
t397 = t283 * t368;
t282 = -sin(t308) + sin(t309);
t396 = t282 * t370;
t336 = sin(pkin(7));
t395 = pkin(1) * t336;
t394 = Ifges(3,3) * t351;
t393 = m(3) * pkin(2) + mrSges(2,1);
t323 = qJ(1,3) + t338;
t324 = qJ(1,3) - t338;
t267 = t282 * pkin(3) + (-sin(t314) + sin(t315)) * pkin(2) + (-sin(t323) + sin(t324)) * pkin(1);
t392 = t267 * t290;
t325 = qJ(1,2) + t339;
t326 = qJ(1,2) - t339;
t268 = t283 * pkin(3) + (-sin(t316) + sin(t317)) * pkin(2) + (-sin(t325) + sin(t326)) * pkin(1);
t391 = t268 * t291;
t327 = qJ(1,1) + t340;
t328 = qJ(1,1) - t340;
t269 = t284 * pkin(3) + (-sin(t318) + sin(t319)) * pkin(2) + (-sin(t327) + sin(t328)) * pkin(1);
t390 = t269 * t292;
t270 = -t285 * pkin(3) + (-cos(t315) - cos(t314)) * pkin(2) + (-cos(t324) - cos(t323)) * pkin(1);
t389 = t270 * t290;
t271 = -t286 * pkin(3) + (-cos(t317) - cos(t316)) * pkin(2) + (-cos(t326) - cos(t325)) * pkin(1);
t388 = t271 * t291;
t272 = -t287 * pkin(3) + (-cos(t319) - cos(t318)) * pkin(2) + (-cos(t328) - cos(t327)) * pkin(1);
t387 = t272 * t292;
t337 = cos(pkin(7));
t356 = pkin(1) * t337 + pkin(2);
t288 = -mrSges(3,1) * t395 - t356 * mrSges(3,2);
t289 = t356 * mrSges(3,1) - mrSges(3,2) * t395;
t344 = cos(qJ(3,3));
t276 = t288 * t341 + t289 * t344 + Ifges(3,3);
t386 = t276 * t351;
t345 = cos(qJ(3,2));
t277 = t288 * t342 + t289 * t345 + Ifges(3,3);
t385 = t277 * t351;
t346 = cos(qJ(3,1));
t278 = t288 * t343 + t289 * t346 + Ifges(3,3);
t384 = t278 * t351;
t320 = sin(qJ(1,3) + t358);
t279 = pkin(3) * t320 + pkin(2) * sin(t361) + sin(qJ(1,3)) * pkin(1);
t383 = t279 * t290;
t382 = t279 * t351;
t321 = sin(qJ(1,2) + t359);
t280 = pkin(3) * t321 + pkin(2) * sin(t362) + sin(qJ(1,2)) * pkin(1);
t381 = t280 * t291;
t380 = t280 * t351;
t322 = sin(qJ(1,1) + t360);
t281 = pkin(3) * t322 + pkin(2) * sin(t363) + sin(qJ(1,1)) * pkin(1);
t379 = t281 * t292;
t378 = t281 * t351;
t371 = t290 * t320;
t369 = t291 * t321;
t367 = t292 * t322;
t365 = 0.2e1 * pkin(1);
t364 = 0.2e1 * pkin(2);
t329 = sin(t338);
t330 = sin(t339);
t331 = sin(t340);
t332 = cos(t338);
t333 = cos(t339);
t334 = cos(t340);
t348 = m(2) + m(3);
t357 = (t329 * t332 + t330 * t333 + t331 * t334) * t348;
t355 = m(3) * pkin(2) ^ 2 + t348 * pkin(1) ^ 2 + Ifges(1,3) + Ifges(2,3) + Ifges(3,3);
t354 = t344 * mrSges(3,1) - mrSges(3,2) * t341;
t353 = t345 * mrSges(3,1) - mrSges(3,2) * t342;
t352 = t346 * mrSges(3,1) - mrSges(3,2) * t343;
t275 = t352 * t364 + ((t352 + t393) * t337 - (mrSges(3,1) * t343 + t346 * mrSges(3,2) + mrSges(2,2)) * t336) * t365 + t355;
t274 = t353 * t364 + ((t353 + t393) * t337 - (mrSges(3,1) * t342 + t345 * mrSges(3,2) + mrSges(2,2)) * t336) * t365 + t355;
t273 = t354 * t364 + ((t354 + t393) * t337 - (mrSges(3,1) * t341 + t344 * mrSges(3,2) + mrSges(2,2)) * t336) * t365 + t355;
t266 = (Ifges(3,3) * t378 - t278 * t322) * t292;
t265 = (Ifges(3,3) * t380 - t277 * t321) * t291;
t264 = (Ifges(3,3) * t382 - t276 * t320) * t290;
t263 = (-t275 * t322 + t278 * t378) * t292;
t262 = (-t274 * t321 + t277 * t380) * t291;
t261 = (-t273 * t320 + t276 * t382) * t290;
t260 = (t272 * t394 + t278 * t287) * t366;
t259 = (t271 * t394 + t277 * t286) * t368;
t258 = (t270 * t394 + t276 * t285) * t370;
t257 = (-t269 * t394 + t278 * t284) * t366;
t256 = (-t268 * t394 + t277 * t283) * t368;
t255 = (-t267 * t394 + t276 * t282) * t370;
t254 = (t272 * t384 + t275 * t287) * t366;
t253 = (t271 * t385 + t274 * t286) * t368;
t252 = (t270 * t386 + t273 * t285) * t370;
t251 = (-t269 * t384 + t275 * t284) * t366;
t250 = (-t268 * t385 + t274 * t283) * t368;
t249 = (-t267 * t386 + t273 * t282) * t370;
t1 = [m(4) + (t329 ^ 2 + t330 ^ 2 + t331 ^ 2) * t348 + t252 * t399 + t253 * t400 + t254 * t401 + (t258 * t389 + t259 * t388 + t260 * t387) * t402, t252 * t396 + t253 * t397 + t254 * t398 + (-t258 * t392 - t259 * t391 - t260 * t390) * t402 + t357, -t252 * t371 - t253 * t369 - t254 * t367 + (t258 * t383 + t259 * t381 + t260 * t379) * t351; t249 * t399 + t250 * t400 + t251 * t401 + (t255 * t389 + t256 * t388 + t257 * t387) * t402 + t357, m(4) + (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) * t348 + t249 * t396 + t250 * t397 + t251 * t398 + (-t255 * t392 - t256 * t391 - t257 * t390) * t402, -t249 * t371 - t250 * t369 - t251 * t367 + (t255 * t383 + t256 * t381 + t257 * t379) * t351; t261 * t399 + t262 * t400 + t263 * t401 + (t264 * t389 + t265 * t388 + t266 * t387) * t402, t261 * t396 + t262 * t397 + t263 * t398 + (-t264 * t392 - t265 * t391 - t266 * t390) * t402, -t261 * t371 - t262 * t369 - t263 * t367 + m(4) + (t264 * t383 + t265 * t381 + t266 * t379) * t351;];
MX  = t1;
