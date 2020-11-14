% Calculate inertia matrix for parallel robot
% P3PRRRR1G2A0
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
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:28
% EndTime: 2020-03-09 21:16:28
% DurationCPUTime: 0.48s
% Computational Cost: add. (702->124), mult. (1140->244), div. (576->10), fcn. (1497->18), ass. (0->124)
t369 = 2 * Ifges(3,4);
t368 = Ifges(3,1) + Ifges(2,3);
t301 = sin(qJ(3,3));
t307 = cos(qJ(3,3));
t277 = mrSges(3,1) * t301 + mrSges(3,2) * t307;
t290 = 0.1e1 / t307;
t367 = t277 * t290;
t303 = sin(qJ(3,2));
t309 = cos(qJ(3,2));
t278 = mrSges(3,1) * t303 + mrSges(3,2) * t309;
t292 = 0.1e1 / t309;
t366 = t278 * t292;
t305 = sin(qJ(3,1));
t311 = cos(qJ(3,1));
t279 = mrSges(3,1) * t305 + mrSges(3,2) * t311;
t294 = 0.1e1 / t311;
t365 = t279 * t294;
t298 = legFrame(3,2);
t280 = sin(t298);
t313 = 0.1e1 / pkin(2);
t364 = t280 * t313;
t299 = legFrame(2,2);
t281 = sin(t299);
t363 = t281 * t313;
t300 = legFrame(1,2);
t282 = sin(t300);
t362 = t282 * t313;
t283 = cos(t298);
t361 = t283 * t313;
t284 = cos(t299);
t360 = t284 * t313;
t285 = cos(t300);
t359 = t285 * t313;
t302 = sin(qJ(2,3));
t287 = 0.1e1 / t302;
t358 = t287 * t290;
t308 = cos(qJ(2,3));
t357 = t287 * t308;
t304 = sin(qJ(2,2));
t288 = 0.1e1 / t304;
t356 = t288 * t292;
t310 = cos(qJ(2,2));
t355 = t288 * t310;
t306 = sin(qJ(2,1));
t289 = 0.1e1 / t306;
t354 = t289 * t294;
t312 = cos(qJ(2,1));
t353 = t289 * t312;
t352 = t290 * t313;
t351 = t292 * t313;
t350 = t294 * t313;
t349 = t302 * t307;
t348 = t304 * t309;
t347 = t306 * t311;
t265 = t280 * t349 - t283 * t301;
t346 = t265 * t358;
t266 = t281 * t348 - t284 * t303;
t345 = t266 * t356;
t267 = t282 * t347 - t285 * t305;
t344 = t267 * t354;
t268 = t280 * t301 + t283 * t349;
t343 = t268 * t358;
t269 = t281 * t303 + t284 * t348;
t342 = t269 * t356;
t270 = t282 * t305 + t285 * t347;
t341 = t270 * t354;
t286 = m(1) + m(2) + m(3);
t340 = t286 * t358;
t339 = t286 * t356;
t338 = t286 * t354;
t314 = t307 ^ 2;
t337 = 0.1e1 / t314 * t301 * t357;
t315 = t309 ^ 2;
t336 = 0.1e1 / t315 * t303 * t355;
t316 = t311 ^ 2;
t335 = 0.1e1 / t316 * t305 * t353;
t274 = Ifges(3,5) * t301 + Ifges(3,6) * t307;
t334 = -Ifges(3,3) * t290 + t274 * t337;
t275 = Ifges(3,5) * t303 + Ifges(3,6) * t309;
t333 = -Ifges(3,3) * t292 + t275 * t336;
t276 = Ifges(3,5) * t305 + Ifges(3,6) * t311;
t332 = -Ifges(3,3) * t294 + t276 * t335;
t297 = mrSges(2,2) - mrSges(3,3);
t262 = (t307 * mrSges(3,1) - t301 * mrSges(3,2) + mrSges(2,1)) * t308 - t302 * t297;
t296 = -Ifges(3,1) + Ifges(3,2);
t271 = t301 * t307 * t369 + t296 * t314 + t368;
t322 = t271 * t337 - t274 * t290;
t235 = t262 * t346 + t322 * t361;
t331 = t235 * t337 - (-t265 * t367 + t334 * t361) * t290;
t263 = (t309 * mrSges(3,1) - t303 * mrSges(3,2) + mrSges(2,1)) * t310 - t304 * t297;
t272 = t303 * t309 * t369 + t296 * t315 + t368;
t321 = t272 * t336 - t275 * t292;
t236 = t263 * t345 + t321 * t360;
t330 = t236 * t336 - (-t266 * t366 + t333 * t360) * t292;
t264 = (t311 * mrSges(3,1) - t305 * mrSges(3,2) + mrSges(2,1)) * t312 - t306 * t297;
t273 = t305 * t311 * t369 + t296 * t316 + t368;
t320 = t273 * t335 - t276 * t294;
t237 = t264 * t344 + t320 * t359;
t329 = t237 * t335 - (-t267 * t365 + t332 * t359) * t294;
t238 = t262 * t343 - t322 * t364;
t328 = t238 * t337 - (-t268 * t367 - t334 * t364) * t290;
t239 = t263 * t342 - t321 * t363;
t327 = t239 * t336 - (-t269 * t366 - t333 * t363) * t292;
t240 = t264 * t341 - t320 * t362;
t326 = t240 * t335 - (-t270 * t365 - t332 * t362) * t294;
t253 = (t262 * t308 - t271 * t352) * t287;
t325 = t253 * t337 - (-t287 * t274 * t352 - t308 * t277) * t290;
t254 = (t263 * t310 - t272 * t351) * t288;
t324 = t254 * t336 - (-t288 * t275 * t351 - t310 * t278) * t292;
t255 = (t264 * t312 - t273 * t350) * t289;
t323 = t255 * t335 - (-t289 * t276 * t350 - t312 * t279) * t294;
t319 = t262 * t337 + t302 * t367;
t318 = t263 * t336 + t304 * t366;
t317 = t264 * t335 + t306 * t365;
t258 = (-t264 * t350 + t286 * t312) * t289;
t257 = (-t263 * t351 + t286 * t310) * t288;
t256 = (-t262 * t352 + t286 * t308) * t287;
t246 = t270 * t338 - t317 * t362;
t245 = t269 * t339 - t318 * t363;
t244 = t268 * t340 - t319 * t364;
t243 = t267 * t338 + t317 * t359;
t242 = t266 * t339 + t318 * t360;
t241 = t265 * t340 + t319 * t361;
t1 = [t241 * t346 + t242 * t345 + t243 * t344 + m(4) + (t331 * t283 + t330 * t284 + t329 * t285) * t313, t241 * t343 + t242 * t342 + t243 * t341 + (-t331 * t280 - t330 * t281 - t329 * t282) * t313, t241 * t357 + t242 * t355 + t243 * t353 + (-t235 * t358 - t236 * t356 - t237 * t354) * t313; t244 * t346 + t245 * t345 + t246 * t344 + (t328 * t283 + t327 * t284 + t326 * t285) * t313, t244 * t343 + t245 * t342 + t246 * t341 + m(4) + (-t328 * t280 - t327 * t281 - t326 * t282) * t313, t244 * t357 + t245 * t355 + t246 * t353 + (-t238 * t358 - t239 * t356 - t240 * t354) * t313; t256 * t346 + t257 * t345 + t258 * t344 + (t325 * t283 + t324 * t284 + t323 * t285) * t313, t256 * t343 + t257 * t342 + t258 * t341 + (-t325 * t280 - t324 * t281 - t323 * t282) * t313, t256 * t357 + t257 * t355 + t258 * t353 + m(4) + (-t253 * t358 - t254 * t356 - t255 * t354) * t313;];
MX  = t1;
