% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x8]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRR1A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRR1A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1A0_coriolisvec_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1A0_coriolisvec_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:40
% EndTime: 2019-05-03 14:47:40
% DurationCPUTime: 0.58s
% Computational Cost: add. (855->99), mult. (2096->194), div. (369->8), fcn. (1612->14), ass. (0->91)
t305 = xP(3);
t281 = sin(t305);
t282 = cos(t305);
t306 = koppelP(3,2);
t309 = koppelP(3,1);
t269 = t281 * t309 + t282 * t306;
t293 = legFrame(3,3);
t275 = sin(t293);
t278 = cos(t293);
t325 = t281 * t306 - t282 * t309;
t261 = t269 * t278 + t325 * t275;
t307 = koppelP(2,2);
t310 = koppelP(2,1);
t270 = t281 * t310 + t282 * t307;
t294 = legFrame(2,3);
t276 = sin(t294);
t279 = cos(t294);
t324 = t281 * t307 - t282 * t310;
t262 = t270 * t279 + t324 * t276;
t308 = koppelP(1,2);
t311 = koppelP(1,1);
t271 = t281 * t311 + t282 * t308;
t295 = legFrame(1,3);
t277 = sin(t295);
t280 = cos(t295);
t323 = t281 * t308 - t282 * t311;
t260 = t271 * t280 + t323 * t277;
t302 = xDP(3);
t303 = xDP(2);
t304 = xDP(1);
t257 = t261 * t302 - t275 * t303 - t304 * t278;
t296 = sin(qJ(2,3));
t349 = t257 ^ 2 / t296 ^ 2;
t258 = t262 * t302 - t276 * t303 - t304 * t279;
t297 = sin(qJ(2,2));
t347 = t258 ^ 2 / t297 ^ 2;
t259 = t260 * t302 - t277 * t303 - t304 * t280;
t298 = sin(qJ(2,1));
t345 = t259 ^ 2 / t298 ^ 2;
t283 = 0.1e1 / t296;
t286 = 0.1e1 / t297;
t289 = 0.1e1 / t298;
t299 = cos(qJ(2,3));
t263 = -t275 * t296 + t278 * t299;
t264 = t275 * t299 + t296 * t278;
t292 = t302 ^ 2;
t312 = 0.1e1 / pkin(2);
t348 = t283 * t349;
t245 = t312 * t348 + (t263 * t325 - t264 * t269) * t292 * t283;
t355 = t245 * t283;
t300 = cos(qJ(2,2));
t265 = -t276 * t297 + t279 * t300;
t266 = t276 * t300 + t297 * t279;
t346 = t286 * t347;
t246 = t312 * t346 + (t265 * t324 - t266 * t270) * t292 * t286;
t354 = t246 * t286;
t301 = cos(qJ(2,1));
t267 = -t277 * t298 + t280 * t301;
t268 = t277 * t301 + t298 * t280;
t344 = t289 * t345;
t247 = t312 * t344 + (t267 * t323 - t268 * t271) * t292 * t289;
t353 = t247 * t289;
t313 = 0.1e1 / pkin(2) ^ 2;
t328 = t275 * t269 - t278 * t325;
t331 = t299 * t348;
t338 = t292 * t312;
t248 = t328 * t283 * t338 - t313 * t331;
t352 = t248 * t283;
t327 = t276 * t270 - t279 * t324;
t330 = t300 * t346;
t249 = t327 * t286 * t338 - t313 * t330;
t351 = t249 * t286;
t326 = t277 * t271 - t280 * t323;
t329 = t301 * t344;
t250 = t326 * t289 * t338 - t313 * t329;
t350 = t250 * t289;
t343 = t281 * t292;
t342 = t282 * t292;
t341 = t283 * t299;
t340 = t286 * t300;
t339 = t289 * t301;
t337 = t245 * t341;
t336 = t246 * t340;
t335 = t247 * t339;
t334 = t248 * t341;
t333 = t249 * t340;
t332 = t250 * t339;
t253 = -t260 * t301 + t326 * t298;
t252 = -t262 * t300 + t327 * t297;
t251 = -t261 * t299 + t328 * t296;
t1 = [t263 * t355 + t265 * t354 + t267 * t353, (-t278 * t352 - t279 * t351 - t280 * t350) * t312, t263 * t334 + t265 * t333 + t267 * t332 + (-t263 * t349 - t265 * t347 - t267 * t345) * t313 + (-t278 * t337 - t279 * t336 - t280 * t335) * t312, -t263 * t248 - t265 * t249 - t267 * t250 + (-t263 * t331 - t265 * t330 - t267 * t329) * t313 + (t245 * t278 + t246 * t279 + t247 * t280) * t312, 0, -t342, t343, 0; t264 * t355 + t266 * t354 + t268 * t353, (-t275 * t352 - t276 * t351 - t277 * t350) * t312, t264 * t334 + t266 * t333 + t268 * t332 + (-t264 * t349 - t266 * t347 - t268 * t345) * t313 + (-t275 * t337 - t276 * t336 - t277 * t335) * t312, -t264 * t248 - t266 * t249 - t268 * t250 + (-t264 * t331 - t266 * t330 - t268 * t329) * t313 + (t245 * t275 + t246 * t276 + t247 * t277) * t312, 0, -t343, -t342, 0; t251 * t355 + t252 * t354 + t253 * t353, (t260 * t350 + t261 * t352 + t262 * t351) * t312, t251 * t334 + t252 * t333 + t253 * t332 + (-t251 * t349 - t252 * t347 - t253 * t345) * t313 + (t260 * t335 + t261 * t337 + t262 * t336) * t312, -t251 * t248 - t252 * t249 - t253 * t250 + (-t251 * t331 - t252 * t330 - t253 * t329) * t313 + (-t245 * t261 - t246 * t262 - t247 * t260) * t312, 0, 0, 0, 0;];
tau_reg  = t1;
