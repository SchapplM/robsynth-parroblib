% Calculate vector of inverse dynamics forces for parallel robot
% P6RRPRRR14V4G7A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [6x1]
%   Generalized platform coordinates
% xDP [6x1]
%   Generalized platform velocities
% xDDP [6x1]
%   Generalized platform accelerations
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,alpha3,alpha4,d1,d4,theta3]';
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
% tauA [6x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in actuated joint coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2022-11-03 12:07
% Revision: 5de314b2d97380370c92dd342c670b987640c890 (2022-02-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauA = P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: xP has to be [6x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [6 1]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: xDP has to be [6x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [6 1]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: xDDP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: qJ has to be [3x6] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V4G7A0_invdyn_para_qa_slagn_vp2: Koppelpunkt has to be [6x3] (double)');

%% Function calls and calculation
tauX = P6RRPRRR14V4G7A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges);
Jinv = P6RRPRRR14V4G7A0_Jinv(xP, qJ, pkin, koppelP, legFrame);
tauA  = (Jinv') \ tauX;
tauA  = tauA;
