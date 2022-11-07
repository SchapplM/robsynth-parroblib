% Calculate inertia matrix for parallel robot
% P3RRPRR1V1G2A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4]';
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
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR1V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:08:04
% EndTime: 2022-11-04 17:08:04
% DurationCPUTime: 0.56s
% Computational Cost: add. (1542->154), mult. (1635->302), div. (342->7), fcn. (1191->18), ass. (0->131)
t356 = (pkin(1) * m(3));
t355 = 2 * mrSges(3,3);
t354 = -2 * mrSges(3,2) * pkin(1) + 2 * Ifges(2,4) + 2 * Ifges(3,4);
t353 = (Ifges(2,1) + Ifges(3,1));
t352 = (-Ifges(2,6) - Ifges(3,6));
t280 = (-mrSges(3,1) - t356);
t295 = sin(qJ(2,3));
t301 = cos(qJ(2,3));
t312 = -pkin(1) * mrSges(3,3) + Ifges(2,5) + Ifges(3,5);
t256 = (t280 * qJ(3,3) + t312) * t295 - (mrSges(3,2) * qJ(3,3) + t352) * t301;
t289 = pkin(3) + qJ(3,3);
t281 = 1 / t289;
t351 = t256 * t281;
t307 = pkin(2) + pkin(1);
t288 = 1 / t307;
t350 = t256 * t288;
t297 = sin(qJ(2,2));
t303 = cos(qJ(2,2));
t257 = (t280 * qJ(3,2) + t312) * t297 - (mrSges(3,2) * qJ(3,2) + t352) * t303;
t290 = pkin(3) + qJ(3,2);
t282 = 1 / t290;
t349 = t257 * t282;
t348 = t257 * t288;
t299 = sin(qJ(2,1));
t305 = cos(qJ(2,1));
t258 = (t280 * qJ(3,1) + t312) * t299 - (mrSges(3,2) * qJ(3,1) + t352) * t305;
t291 = pkin(3) + qJ(3,1);
t283 = 1 / t291;
t347 = t258 * t283;
t346 = t258 * t288;
t292 = legFrame(3,2);
t274 = sin(t292);
t277 = cos(t292);
t296 = sin(qJ(1,3));
t328 = t296 * t301;
t259 = -t274 * t328 + t295 * t277;
t285 = 0.1e1 / t301;
t345 = t259 * t285;
t293 = legFrame(2,2);
t275 = sin(t293);
t278 = cos(t293);
t298 = sin(qJ(1,2));
t326 = t298 * t303;
t260 = -t275 * t326 + t297 * t278;
t286 = 0.1e1 / t303;
t344 = t260 * t286;
t294 = legFrame(1,2);
t276 = sin(t294);
t279 = cos(t294);
t300 = sin(qJ(1,1));
t324 = t300 * t305;
t261 = -t276 * t324 + t299 * t279;
t287 = 0.1e1 / t305;
t343 = t261 * t287;
t262 = t274 * t295 + t277 * t328;
t342 = t262 * t285;
t263 = t275 * t297 + t278 * t326;
t341 = t263 * t286;
t264 = t276 * t299 + t279 * t324;
t340 = t264 * t287;
t268 = t295 * mrSges(3,2) + t280 * t301;
t339 = t268 * t285;
t269 = t297 * mrSges(3,2) + t280 * t303;
t338 = t269 * t286;
t270 = t299 * mrSges(3,2) + t280 * t305;
t337 = t270 * t287;
t320 = (2 * mrSges(3,1) + t356) * pkin(1);
t336 = (Ifges(2,3) + Ifges(3,3) + t320) * t288;
t335 = t274 * t285;
t334 = t275 * t286;
t333 = t276 * t287;
t332 = t277 * t285;
t331 = t278 * t286;
t330 = t279 * t287;
t329 = t295 * t307;
t327 = t297 * t307;
t325 = t299 * t307;
t323 = t301 * t307;
t322 = t303 * t307;
t321 = t305 * t307;
t319 = Ifges(1,3) + t353;
t318 = t285 * t350;
t302 = cos(qJ(1,3));
t317 = t302 * t350;
t316 = t286 * t348;
t304 = cos(qJ(1,2));
t315 = t304 * t348;
t314 = t287 * t346;
t306 = cos(qJ(1,1));
t313 = t306 * t346;
t311 = -t289 * t302 + t296 * t323;
t310 = -t290 * t304 + t298 * t322;
t309 = -t291 * t306 + t300 * t321;
t271 = Ifges(2,2) + Ifges(3,2) + t320 - t353;
t267 = t300 * t291 + t306 * t321;
t266 = t298 * t290 + t304 * t322;
t265 = t296 * t289 + t302 * t323;
t255 = t276 * t325 + t309 * t279;
t254 = t275 * t327 + t310 * t278;
t253 = t274 * t329 + t311 * t277;
t252 = -t309 * t276 + t279 * t325;
t251 = -t310 * t275 + t278 * t327;
t250 = -t311 * t274 + t277 * t329;
t249 = (t271 * t305 + t299 * t354) * t305 + ((m(3) * qJ(3,1) + t355) * qJ(3,1)) + t319;
t248 = (t271 * t303 + t297 * t354) * t303 + ((m(3) * qJ(3,2) + t355) * qJ(3,2)) + t319;
t247 = (t271 * t301 + t295 * t354) * t301 + ((m(3) * qJ(3,3) + t355) * qJ(3,3)) + t319;
t246 = (m(3) * t267 + t270 * t306) * t283;
t245 = (m(3) * t266 + t269 * t304) * t282;
t244 = (m(3) * t265 + t268 * t302) * t281;
t243 = (t264 * t347 + t276 * t336) * t287;
t242 = (t263 * t349 + t275 * t336) * t286;
t241 = (t262 * t351 + t274 * t336) * t285;
t240 = (t261 * t347 + t279 * t336) * t287;
t239 = (t260 * t349 + t278 * t336) * t286;
t238 = (t259 * t351 + t277 * t336) * t285;
t237 = (t249 * t306 + t267 * t270) * t283;
t236 = (t248 * t304 + t266 * t269) * t282;
t235 = (t247 * t302 + t265 * t268) * t281;
t234 = (m(3) * t255 + t264 * t337) * t283;
t233 = (m(3) * t254 + t263 * t338) * t282;
t232 = (m(3) * t253 + t262 * t339) * t281;
t231 = (m(3) * t252 + t261 * t337) * t283;
t230 = (m(3) * t251 + t260 * t338) * t282;
t229 = (m(3) * t250 + t259 * t339) * t281;
t228 = t276 * t314 + (t249 * t340 + t255 * t270) * t283;
t227 = t275 * t316 + (t248 * t341 + t254 * t269) * t282;
t226 = t274 * t318 + (t247 * t342 + t253 * t268) * t281;
t225 = t279 * t314 + (t249 * t343 + t252 * t270) * t283;
t224 = t278 * t316 + (t248 * t344 + t251 * t269) * t282;
t223 = t277 * t318 + (t247 * t345 + t250 * t268) * t281;
t1 = [m(4) + (t228 * t340 + t234 * t255) * t283 + (t227 * t341 + t233 * t254) * t282 + (t226 * t342 + t232 * t253) * t281 + (t241 * t335 + t242 * t334 + t243 * t333) * t288, (t228 * t343 + t234 * t252) * t283 + (t227 * t344 + t233 * t251) * t282 + (t226 * t345 + t232 * t250) * t281 + (t241 * t332 + t242 * t331 + t243 * t330) * t288, (t228 * t306 + t234 * t267) * t283 + (t227 * t304 + t233 * t266) * t282 + (t226 * t302 + t232 * t265) * t281; (t225 * t340 + t231 * t255) * t283 + (t224 * t341 + t230 * t254) * t282 + (t223 * t342 + t229 * t253) * t281 + (t238 * t335 + t239 * t334 + t240 * t333) * t288, m(4) + (t225 * t343 + t231 * t252) * t283 + (t224 * t344 + t230 * t251) * t282 + (t223 * t345 + t229 * t250) * t281 + (t238 * t332 + t239 * t331 + t240 * t330) * t288, (t225 * t306 + t231 * t267) * t283 + (t224 * t304 + t230 * t266) * t282 + (t223 * t302 + t229 * t265) * t281; (t246 * t255 + (t237 * t264 + t276 * t313) * t287) * t283 + (t245 * t254 + (t236 * t263 + t275 * t315) * t286) * t282 + (t244 * t253 + (t235 * t262 + t274 * t317) * t285) * t281, (t246 * t252 + (t237 * t261 + t279 * t313) * t287) * t283 + (t245 * t251 + (t236 * t260 + t278 * t315) * t286) * t282 + (t244 * t250 + (t235 * t259 + t277 * t317) * t285) * t281, m(4) + (t237 * t306 + t246 * t267) * t283 + (t236 * t304 + t245 * t266) * t282 + (t235 * t302 + t244 * t265) * t281;];
MX  = t1;
