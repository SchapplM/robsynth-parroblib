% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PPR1A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [6x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PPR1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PPR1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PPR1A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1A0_coriolisvec_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1A0_coriolisvec_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'P3PPR1A0_coriolisvec_para_pf_mdp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:32
% EndTime: 2019-05-03 14:37:32
% DurationCPUTime: 0.14s
% Computational Cost: add. (136->46), mult. (369->85), div. (0->0), fcn. (238->8), ass. (0->41)
t371 = xP(3);
t365 = sin(t371);
t366 = cos(t371);
t372 = koppelP(3,2);
t375 = koppelP(3,1);
t353 = t365 * t375 + t366 * t372;
t356 = -t365 * t372 + t366 * t375;
t368 = legFrame(3,3);
t359 = sin(t368);
t362 = cos(t368);
t367 = xDP(3) ^ 2;
t341 = (-t353 * t362 + t356 * t359) * t367;
t373 = koppelP(2,2);
t376 = koppelP(2,1);
t354 = t365 * t376 + t366 * t373;
t357 = -t365 * t373 + t366 * t376;
t369 = legFrame(2,3);
t360 = sin(t369);
t363 = cos(t369);
t342 = (-t354 * t363 + t357 * t360) * t367;
t374 = koppelP(1,2);
t377 = koppelP(1,1);
t355 = t365 * t377 + t366 * t374;
t358 = -t365 * t374 + t366 * t377;
t370 = legFrame(1,3);
t361 = sin(t370);
t364 = cos(t370);
t343 = (-t355 * t364 + t358 * t361) * t367;
t347 = t359 * t375 - t362 * t372;
t348 = t359 * t372 + t362 * t375;
t349 = t360 * t376 - t363 * t373;
t350 = t360 * t373 + t363 * t376;
t351 = t361 * t377 - t364 * t374;
t352 = t361 * t374 + t364 * t377;
t380 = (t365 * t351 + t352 * t366) * t343 + (t365 * t349 + t350 * t366) * t342 + (t365 * t347 + t348 * t366) * t341;
t379 = t362 * t341 + t363 * t342 + t364 * t343;
t378 = -t359 * t341 - t360 * t342 - t361 * t343;
t346 = (-t355 * t361 - t358 * t364) * t367;
t345 = (-t354 * t360 - t357 * t363) * t367;
t344 = (-t353 * t359 - t356 * t362) * t367;
t1 = [t378 * MDP(1) + (t362 * t344 + t363 * t345 + t364 * t346 + t378) * MDP(2) + (-t366 * MDP(4) + t365 * MDP(5)) * t367; t379 * MDP(1) + (t359 * t344 + t360 * t345 + t361 * t346 + t379) * MDP(2) + (-t365 * MDP(4) - t366 * MDP(5)) * t367; t380 * MDP(1) + ((t351 * t366 - t352 * t365) * t346 + (t349 * t366 - t350 * t365) * t345 + (t347 * t366 - t348 * t365) * t344 + t380) * MDP(2);];
taucX  = t1;
