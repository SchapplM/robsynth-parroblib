% Calculate inertia matrix for parallel robot
% P3PRRRR1G1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 20:34
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 20:34:06
% EndTime: 2020-03-09 20:34:06
% DurationCPUTime: 0.43s
% Computational Cost: add. (561->116), mult. (951->235), div. (408->10), fcn. (1272->18), ass. (0->145)
t404 = 2 * Ifges(3,4);
t403 = Ifges(3,1) + Ifges(2,3);
t335 = 0.1e1 / pkin(2);
t402 = Ifges(3,3) * t335;
t322 = mrSges(2,2) - mrSges(3,3);
t323 = sin(qJ(3,3));
t324 = sin(qJ(2,3));
t329 = cos(qJ(3,3));
t330 = cos(qJ(2,3));
t275 = (t329 * mrSges(3,1) - t323 * mrSges(3,2) + mrSges(2,1)) * t330 - t324 * t322;
t312 = 0.1e1 / t329;
t401 = t275 * t312;
t325 = sin(qJ(3,2));
t326 = sin(qJ(2,2));
t331 = cos(qJ(3,2));
t332 = cos(qJ(2,2));
t276 = (t331 * mrSges(3,1) - t325 * mrSges(3,2) + mrSges(2,1)) * t332 - t326 * t322;
t314 = 0.1e1 / t331;
t400 = t276 * t314;
t327 = sin(qJ(3,1));
t328 = sin(qJ(2,1));
t333 = cos(qJ(3,1));
t334 = cos(qJ(2,1));
t277 = (t333 * mrSges(3,1) - t327 * mrSges(3,2) + mrSges(2,1)) * t334 - t328 * t322;
t316 = 0.1e1 / t333;
t399 = t277 * t316;
t318 = legFrame(3,3);
t302 = sin(t318);
t398 = t302 * t312;
t319 = legFrame(2,3);
t303 = sin(t319);
t397 = t303 * t314;
t320 = legFrame(1,3);
t304 = sin(t320);
t396 = t304 * t316;
t305 = cos(t318);
t395 = t305 * t312;
t306 = cos(t319);
t394 = t306 * t314;
t307 = cos(t320);
t393 = t307 * t316;
t309 = 0.1e1 / t324;
t392 = t309 * t312;
t336 = t329 ^ 2;
t313 = 0.1e1 / t336;
t391 = t309 * t313;
t310 = 0.1e1 / t326;
t390 = t310 * t314;
t337 = t331 ^ 2;
t315 = 0.1e1 / t337;
t389 = t310 * t315;
t311 = 0.1e1 / t328;
t388 = t311 * t316;
t338 = t333 ^ 2;
t317 = 0.1e1 / t338;
t387 = t311 * t317;
t386 = t313 * t335;
t385 = t315 * t335;
t384 = t317 * t335;
t383 = t323 * t330;
t382 = t325 * t332;
t381 = t327 * t334;
t380 = t329 * t330;
t379 = t331 * t332;
t378 = t333 * t334;
t281 = -t302 * t383 - t305 * t329;
t284 = t302 * t323 + t305 * t380;
t350 = t309 * t386;
t347 = t275 * t350;
t308 = m(1) + m(2) + m(3);
t353 = t308 * t392;
t377 = t281 * t347 + t284 * t353;
t282 = -t303 * t382 - t306 * t331;
t286 = t303 * t325 + t306 * t379;
t349 = t310 * t385;
t346 = t276 * t349;
t352 = t308 * t390;
t376 = t282 * t346 + t286 * t352;
t283 = -t304 * t381 - t307 * t333;
t288 = t304 * t327 + t307 * t378;
t348 = t311 * t384;
t345 = t277 * t348;
t351 = t308 * t388;
t375 = t283 * t345 + t288 * t351;
t374 = t281 * t391;
t373 = t282 * t389;
t372 = t283 * t387;
t371 = t284 * t392;
t285 = -t329 * t302 + t305 * t383;
t370 = t285 * t391;
t369 = t286 * t390;
t287 = -t331 * t303 + t306 * t382;
t368 = t287 * t389;
t367 = t288 * t388;
t289 = -t333 * t304 + t307 * t381;
t366 = t289 * t387;
t290 = t302 * t380 - t305 * t323;
t365 = t290 * t392;
t291 = t303 * t379 - t306 * t325;
t364 = t291 * t390;
t292 = t304 * t378 - t307 * t327;
t363 = t292 * t388;
t321 = -Ifges(3,1) + Ifges(3,2);
t362 = (t323 * t329 * t404 + t321 * t336 + t403) * t386;
t361 = (t325 * t331 * t404 + t321 * t337 + t403) * t385;
t360 = (t327 * t333 * t404 + t321 * t338 + t403) * t384;
t296 = Ifges(3,5) * t323 + Ifges(3,6) * t329;
t359 = t296 * t312 * t335;
t297 = Ifges(3,5) * t325 + Ifges(3,6) * t331;
t358 = t297 * t314 * t335;
t298 = Ifges(3,5) * t327 + Ifges(3,6) * t333;
t357 = t298 * t316 * t335;
t299 = mrSges(3,1) * t323 + mrSges(3,2) * t329;
t356 = t299 * t312 * t324;
t300 = mrSges(3,1) * t325 + mrSges(3,2) * t331;
t355 = t300 * t314 * t326;
t301 = mrSges(3,1) * t327 + mrSges(3,2) * t333;
t354 = t301 * t316 * t328;
t341 = t335 * t356;
t254 = t285 * t347 + t290 * t353 + t305 * t341;
t340 = t335 * t355;
t255 = t287 * t346 + t291 * t352 + t306 * t340;
t339 = t335 * t354;
t256 = t289 * t345 + t292 * t351 + t307 * t339;
t344 = t296 * t350;
t343 = t297 * t349;
t342 = t298 * t348;
t262 = t289 * t342 + (-t292 * t301 - t307 * t402) * t316;
t261 = t287 * t343 + (-t291 * t300 - t306 * t402) * t314;
t260 = t285 * t344 + (-t290 * t299 - t305 * t402) * t312;
t259 = t283 * t342 + (-t288 * t301 + t304 * t402) * t316;
t258 = t282 * t343 + (-t286 * t300 + t303 * t402) * t314;
t257 = t281 * t344 + (-t284 * t299 + t302 * t402) * t312;
t253 = -t304 * t339 + t375;
t252 = -t303 * t340 + t376;
t251 = -t302 * t341 + t377;
t250 = -t307 * t357 + (t289 * t360 + t292 * t399) * t311;
t249 = -t306 * t358 + (t287 * t361 + t291 * t400) * t310;
t248 = -t305 * t359 + (t285 * t362 + t290 * t401) * t309;
t247 = t304 * t357 + (t283 * t360 + t288 * t399) * t311;
t246 = t303 * t358 + (t282 * t361 + t286 * t400) * t310;
t245 = t302 * t359 + (t281 * t362 + t284 * t401) * t309;
t244 = t256 + t255 + t254;
t243 = (-t302 * t356 - t303 * t355 - t304 * t354) * t335 + t375 + t376 + t377;
t1 = [t251 * t371 + t252 * t369 + t253 * t367 + m(4) + (t245 * t374 + t246 * t373 + t247 * t372 + t257 * t398 + t258 * t397 + t259 * t396) * t335, t251 * t365 + t252 * t364 + t253 * t363 + (t245 * t370 + t246 * t368 + t247 * t366 - t257 * t395 - t258 * t394 - t259 * t393) * t335, t243; t254 * t371 + t255 * t369 + t256 * t367 + (t248 * t374 + t249 * t373 + t250 * t372 + t260 * t398 + t261 * t397 + t262 * t396) * t335, t254 * t365 + t255 * t364 + t256 * t363 + m(4) + (t248 * t370 + t249 * t368 + t250 * t366 - t260 * t395 - t261 * t394 - t262 * t393) * t335, t244; t243, t244, (3 * m(1)) + (3 * m(2)) + (3 * m(3)) + m(4);];
MX  = t1;
