% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRRR1G3P1A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x12]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR1G3P1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3P1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:25
% EndTime: 2020-03-09 21:02:27
% DurationCPUTime: 1.52s
% Computational Cost: add. (442->172), mult. (1686->399), div. (1002->17), fcn. (2034->18), ass. (0->218)
t289 = legFrame(3,2);
t259 = cos(t289);
t290 = legFrame(2,2);
t260 = cos(t290);
t291 = legFrame(1,2);
t261 = cos(t291);
t303 = cos(qJ(2,1));
t296 = sin(qJ(3,1));
t271 = t296 ^ 2;
t297 = sin(qJ(2,1));
t273 = 0.1e1 / t297 ^ 2;
t302 = cos(qJ(3,1));
t318 = t302 ^ 2;
t285 = 0.1e1 / t318;
t415 = t273 * t285;
t381 = t271 * t415;
t350 = t303 * t381;
t301 = cos(qJ(2,2));
t294 = sin(qJ(3,2));
t267 = t294 ^ 2;
t295 = sin(qJ(2,2));
t269 = 0.1e1 / t295 ^ 2;
t300 = cos(qJ(3,2));
t315 = t300 ^ 2;
t280 = 0.1e1 / t315;
t420 = t269 * t280;
t389 = t267 * t420;
t354 = t301 * t389;
t299 = cos(qJ(2,3));
t292 = sin(qJ(3,3));
t263 = t292 ^ 2;
t293 = sin(qJ(2,3));
t265 = 0.1e1 / t293 ^ 2;
t298 = cos(qJ(3,3));
t312 = t298 ^ 2;
t275 = 0.1e1 / t312;
t425 = t265 * t275;
t397 = t263 * t425;
t358 = t299 * t397;
t448 = t259 * t358 + t260 * t354 + t261 * t350;
t256 = sin(t289);
t447 = t256 * t259;
t257 = sin(t290);
t446 = t257 * t260;
t258 = sin(t291);
t445 = t258 * t261;
t274 = 0.1e1 / t298;
t279 = 0.1e1 / t300;
t284 = 0.1e1 / t302;
t304 = 1 / pkin(2);
t444 = 2 * t304;
t305 = 1 / (pkin(2) ^ 2);
t443 = 2 * t305;
t433 = t256 * t299;
t244 = t259 * t293 - t433;
t264 = 0.1e1 / t293;
t442 = t244 * t264;
t430 = t259 * t299;
t245 = t256 * t293 + t430;
t441 = t245 * t264;
t440 = t245 * t265;
t432 = t257 * t301;
t246 = t260 * t295 - t432;
t268 = 0.1e1 / t295;
t439 = t246 * t268;
t429 = t260 * t301;
t247 = t257 * t295 + t429;
t438 = t247 * t268;
t437 = t247 * t269;
t431 = t258 * t303;
t248 = t261 * t297 - t431;
t272 = 0.1e1 / t297;
t436 = t248 * t272;
t428 = t261 * t303;
t249 = t258 * t297 + t428;
t435 = t249 * t272;
t434 = t249 * t273;
t427 = t263 * t275;
t426 = t264 * t274;
t278 = t299 ^ 2;
t424 = t265 * t278;
t423 = t265 * t299;
t422 = t267 * t280;
t421 = t268 * t279;
t283 = t301 ^ 2;
t419 = t269 * t283;
t418 = t269 * t301;
t417 = t271 * t285;
t416 = t272 * t284;
t288 = t303 ^ 2;
t414 = t273 * t288;
t413 = t273 * t303;
t412 = t274 * t292;
t411 = t275 * t292;
t276 = t274 * t275;
t410 = t276 * t299;
t409 = t279 * t294;
t408 = t280 * t294;
t281 = t279 * t280;
t407 = t281 * t301;
t406 = t284 * t296;
t405 = t285 * t296;
t286 = t284 * t285;
t404 = t286 * t303;
t403 = t256 * t426;
t402 = t257 * t421;
t401 = t258 * t416;
t400 = t259 * t426;
t399 = t260 * t421;
t398 = t261 * t416;
t396 = t278 * t427;
t395 = t264 * t412;
t394 = t264 * t411;
t393 = t265 * t412;
t392 = t276 * t424;
t391 = 0.1e1 / t312 ^ 2 * t424;
t390 = t292 * t423;
t388 = t283 * t422;
t387 = t268 * t409;
t386 = t268 * t408;
t385 = t269 * t409;
t384 = t281 * t419;
t383 = 0.1e1 / t315 ^ 2 * t419;
t382 = t294 * t418;
t380 = t288 * t417;
t379 = t272 * t406;
t378 = t272 * t405;
t377 = t273 * t406;
t376 = t286 * t414;
t375 = 0.1e1 / t318 ^ 2 * t414;
t374 = t296 * t413;
t349 = t284 * t374;
t345 = t258 * t349;
t353 = t279 * t382;
t346 = t257 * t353;
t357 = t274 * t390;
t347 = t256 * t357;
t373 = (t345 + t346 + t347) * t304;
t372 = t448 * t304;
t371 = -0.1e1 - t424;
t370 = -0.1e1 - t419;
t369 = -0.1e1 - t414;
t368 = t244 * t256 * t423;
t367 = t245 * t259 * t423;
t366 = t246 * t257 * t418;
t365 = t247 * t260 * t418;
t364 = t248 * t258 * t413;
t363 = t249 * t261 * t413;
t362 = t425 * t447;
t361 = t420 * t446;
t360 = t415 * t445;
t262 = t292 * t263;
t359 = t262 * t265 * t410;
t356 = t276 * t390;
t266 = t294 * t267;
t355 = t266 * t269 * t407;
t352 = t281 * t382;
t270 = t296 * t271;
t351 = t270 * t273 * t404;
t348 = t286 * t374;
t343 = t259 * t357;
t341 = t260 * t353;
t339 = t261 * t349;
t338 = t244 * t259 - t245 * t256;
t337 = t244 * t278 - t433;
t336 = t246 * t260 - t247 * t257;
t335 = t246 * t283 - t432;
t334 = t248 * t261 - t249 * t258;
t333 = t248 * t288 - t431;
t332 = t284 * t334;
t331 = t338 * t274;
t330 = (-t245 * t278 - t430) * t265;
t329 = t336 * t279;
t328 = (-t247 * t283 - t429) * t269;
t327 = (-t249 * t288 - t428) * t273;
t326 = t338 * t423;
t325 = t336 * t418;
t324 = t334 * t413;
t323 = t262 * t392 + t266 * t384 + t270 * t376;
t322 = t263 * t264 * t410 + t267 * t268 * t407 + t271 * t272 * t404;
t321 = -t256 * t358 - t257 * t354 - t258 * t350;
t255 = t261 ^ 2;
t254 = t260 ^ 2;
t253 = t259 ^ 2;
t252 = t258 ^ 2;
t251 = t257 ^ 2;
t250 = t256 ^ 2;
t237 = (t272 * t380 - t297) * t304;
t236 = (t268 * t388 - t295) * t304;
t235 = (t264 * t396 - t293) * t304;
t234 = (-t272 * t288 - t297) * t304 * t406;
t233 = (-t268 * t283 - t295) * t304 * t409;
t232 = (-t264 * t278 - t293) * t304 * t412;
t231 = (-t398 - t399 - t400) * t305;
t230 = (t401 + t402 + t403) * t305;
t229 = (-t259 * t394 - t260 * t386 - t261 * t378) * t305;
t228 = (t256 * t394 + t257 * t386 + t258 * t378) * t305;
t227 = (-t360 - t361 - t362) * t305;
t226 = (t259 * t356 + t260 * t352 + t261 * t348) * t305;
t225 = (t259 * t359 + t260 * t355 + t261 * t351) * t305;
t224 = (-t256 * t356 - t257 * t352 - t258 * t348) * t305;
t223 = (-t256 * t359 - t257 * t355 - t258 * t351) * t305;
t222 = t448 * t443;
t221 = t321 * t443;
t220 = (-t263 * t362 - t267 * t361 - t271 * t360) * t305;
t219 = (-t377 * t445 - t385 * t446 - t393 * t447) * t443;
t218 = t245 * t393 + t247 * t385 + t249 * t377;
t217 = t244 * t393 + t246 * t385 + t248 * t377;
t216 = t244 * t440 + t246 * t437 + t248 * t434;
t215 = ((t249 * t303 + t261) * t378 + (t247 * t301 + t260) * t386 + (t245 * t299 + t259) * t394) * t304;
t214 = ((t248 * t303 - t258) * t378 + (t246 * t301 - t257) * t386 + (t244 * t299 - t256) * t394) * t304;
t213 = (t327 * t405 + t328 * t408 + t330 * t411) * t304;
t212 = (-t337 * t265 * t411 - t335 * t269 * t408 - t333 * t273 * t405) * t304;
t211 = (t264 * t331 + t268 * t329 + t272 * t332) * t304;
t210 = (-t324 - t325 - t326) * t304;
t209 = (-t274 * t326 - t279 * t325 - t284 * t324) * t304;
t208 = (t329 * t382 + t331 * t390 + t332 * t374) * t304;
t1 = [t245 ^ 2 * t265 + t247 ^ 2 * t269 + t249 ^ 2 * t273, (t253 * t425 + t254 * t420 + t255 * t415) * t305, (-t274 * t367 - t279 * t365 - t284 * t363) * t444, (t245 * t400 + t247 * t399 + t249 * t398) * t444, (t253 * t397 + t254 * t389 + t255 * t381) * t305, (t253 * t393 + t254 * t385 + t255 * t377) * t443, 0, 0, 0, (-t363 - t365 - t367) * t444, (t245 * t343 + t247 * t341 + t249 * t339) * t444, 1; t216, t227, t209, t211, t220, t219, 0, 0, 0, t210, t208, 0; t218, t226, t213, t215, t225, t222, t229, t231, 0, t232 * t441 + t233 * t438 + t234 * t435 + (-t339 - t341 - t343) * t304, t235 * t441 + t236 * t438 + t237 * t435 + t372, 0; t216, t227, t209, t211, t220, t219, 0, 0, 0, t210, t208, 0; t244 ^ 2 * t265 + t246 ^ 2 * t269 + t248 ^ 2 * t273, (t250 * t425 + t251 * t420 + t252 * t415) * t305, (t274 * t368 + t279 * t366 + t284 * t364) * t444, (-t244 * t403 - t246 * t402 - t248 * t401) * t444, (t250 * t397 + t251 * t389 + t252 * t381) * t305, (t250 * t393 + t251 * t385 + t252 * t377) * t443, 0, 0, 0, (t364 + t366 + t368) * t444, (-t244 * t347 - t246 * t346 - t248 * t345) * t444, 1; t217, t224, t212, t214, t223, t221, t228, t230, 0, t232 * t442 + t233 * t439 + t234 * t436 + t373, t235 * t442 + t236 * t439 + t237 * t436 + t304 * t321, 0; t218, t226, t213, t215, t225, t222, t229, t231, 0, ((-t249 + t327) * t406 + (-t247 + t328) * t409 + (-t245 + t330) * t412) * t304, (t380 * t434 + t388 * t437 + t396 * t440 - t245 - t247 - t249) * t304 + t372, 0; t217, t224, t212, t214, t223, t221, t228, t230, 0, (t371 * t244 * t412 + t370 * t246 * t409 + t369 * t248 * t406) * t304 + t373, (t333 * t381 + t335 * t389 + t337 * t397 - t244 - t246 - t248) * t304, 0; t381 + t389 + t397, (t263 * t391 + t267 * t383 + t271 * t375) * t305, (-t263 * t392 - t267 * t384 - t271 * t376) * t444, t322 * t444, (t263 ^ 2 * t391 + t267 ^ 2 * t383 + t271 ^ 2 * t375) * t305, t323 * t443, -0.2e1 * t322 * t305, (-t299 * t394 - t301 * t386 - t303 * t378) * t443, (t275 + t280 + t285) * t305, t232 * t395 + t233 * t387 + t234 * t379 + (t369 * t417 + t370 * t422 + t371 * t427) * t304, t235 * t395 + t236 * t387 + t237 * t379 + (t323 - t406 - t409 - t412) * t304, 1;];
tau_reg  = t1;
