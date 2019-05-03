% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3PPR1A0
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

% Output:
% tau_reg [3x6]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PPR1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1A0_gravload_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PPR1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1A0_gravload_para_pf_regmin: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:32
% EndTime: 2019-05-03 14:37:32
% DurationCPUTime: 0.07s
% Computational Cost: add. (83->37), mult. (151->62), div. (0->0), fcn. (142->8), ass. (0->36)
t199 = legFrame(3,3);
t191 = sin(t199);
t194 = cos(t199);
t203 = koppelP(3,2);
t206 = koppelP(3,1);
t177 = t191 * t206 - t194 * t203;
t178 = t191 * t203 + t194 * t206;
t200 = legFrame(2,3);
t192 = sin(t200);
t195 = cos(t200);
t204 = koppelP(2,2);
t207 = koppelP(2,1);
t179 = t192 * t207 - t195 * t204;
t180 = t192 * t204 + t195 * t207;
t201 = legFrame(1,3);
t193 = sin(t201);
t196 = cos(t201);
t205 = koppelP(1,2);
t208 = koppelP(1,1);
t181 = t193 * t208 - t196 * t205;
t182 = t193 * t205 + t196 * t208;
t183 = t191 * g(1) - t194 * g(2);
t184 = t192 * g(1) - t195 * g(2);
t185 = t193 * g(1) - t196 * g(2);
t202 = xP(3);
t197 = sin(t202);
t198 = cos(t202);
t211 = (t197 * t181 + t182 * t198) * t185 + (t197 * t179 + t180 * t198) * t184 + (t197 * t177 + t178 * t198) * t183;
t210 = t194 * t183 + t195 * t184 + t196 * t185;
t209 = -t191 * t183 - t192 * t184 - t193 * t185;
t190 = t198 * g(1) + t197 * g(2);
t189 = t197 * g(1) - t198 * g(2);
t188 = -t196 * g(1) - t193 * g(2);
t187 = -t195 * g(1) - t192 * g(2);
t186 = -t194 * g(1) - t191 * g(2);
t1 = [t209, t194 * t186 + t195 * t187 + t196 * t188 + t209, 0, 0, 0, -t197 * t189 - t198 * t190; t210, t191 * t186 + t192 * t187 + t193 * t188 + t210, 0, 0, 0, t198 * t189 - t197 * t190; t211, (t181 * t198 - t182 * t197) * t188 + (t179 * t198 - t180 * t197) * t187 + (t177 * t198 - t178 * t197) * t186 + t211, 0, t189, t190, 0;];
tau_reg  = t1;
