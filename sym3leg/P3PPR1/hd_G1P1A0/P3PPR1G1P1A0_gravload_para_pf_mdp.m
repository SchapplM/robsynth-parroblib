% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PPR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
%   see P3PPR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3PPR1G1P1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1G1P1A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1G1P1A0_gravload_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PPR1G1P1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1G1P1A0_gravload_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1G1P1A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1G1P1A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'P3PPR1G1P1A0_gravload_para_pf_mdp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:32
% EndTime: 2019-05-03 14:37:32
% DurationCPUTime: 0.13s
% Computational Cost: add. (90->44), mult. (161->72), div. (0->0), fcn. (142->8), ass. (0->36)
t409 = legFrame(3,3);
t401 = sin(t409);
t404 = cos(t409);
t413 = koppelP(3,2);
t416 = koppelP(3,1);
t387 = t401 * t416 - t404 * t413;
t388 = t401 * t413 + t404 * t416;
t410 = legFrame(2,3);
t402 = sin(t410);
t405 = cos(t410);
t414 = koppelP(2,2);
t417 = koppelP(2,1);
t389 = t402 * t417 - t405 * t414;
t390 = t402 * t414 + t405 * t417;
t411 = legFrame(1,3);
t403 = sin(t411);
t406 = cos(t411);
t415 = koppelP(1,2);
t418 = koppelP(1,1);
t391 = t403 * t418 - t406 * t415;
t392 = t403 * t415 + t406 * t418;
t393 = t401 * g(1) - t404 * g(2);
t394 = t402 * g(1) - t405 * g(2);
t395 = t403 * g(1) - t406 * g(2);
t412 = xP(3);
t407 = sin(t412);
t408 = cos(t412);
t421 = (t407 * t391 + t392 * t408) * t395 + (t407 * t389 + t390 * t408) * t394 + (t407 * t387 + t388 * t408) * t393;
t420 = t404 * t393 + t405 * t394 + t406 * t395;
t419 = -t401 * t393 - t402 * t394 - t403 * t395;
t400 = t408 * g(1) + t407 * g(2);
t399 = t407 * g(1) - t408 * g(2);
t398 = -t406 * g(1) - t403 * g(2);
t397 = -t405 * g(1) - t402 * g(2);
t396 = -t404 * g(1) - t401 * g(2);
t1 = [t419 * MDP(1) + (t404 * t396 + t405 * t397 + t406 * t398 + t419) * MDP(2) + (-t407 * t399 - t408 * t400) * MDP(6); t420 * MDP(1) + (t401 * t396 + t402 * t397 + t403 * t398 + t420) * MDP(2) + (t408 * t399 - t407 * t400) * MDP(6); t421 * MDP(1) + ((t391 * t408 - t392 * t407) * t398 + (t389 * t408 - t390 * t407) * t397 + (t387 * t408 - t388 * t407) * t396 + t421) * MDP(2) + t399 * MDP(4) + t400 * MDP(5);];
taugX  = t1;
